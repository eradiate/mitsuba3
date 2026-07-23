from __future__ import annotations

import drjit as dr
import mitsuba as mi

from .util import SceneParameters, _jit_id_hash


class EradiateSceneParameters(SceneParameters):
    """
    An extended :py:class:`~mitsuba.SceneParameters` that traverses all
    parents of a node. Also collects the nodes' deduplicated ``names``.

    Plain :py:func:`mitsuba.traverse` only remembers the first parent it
    encounters in traversal. :py:meth:`~mitsuba.SceneParameters.update` only
    notifies this parent and ignores any additional parents it might have.
    Additionally, the name of the node it traverses is not guaranteed to be
    based on its id depending on the order it is traversed.

    Use :py:func:`eradiate.traverse` to get an :py:class:`EradiateSceneParameters`
    instead, whose :py:meth:`update` notifies every parent.
    """

    def __init__(self, properties=None, hierarchy=None, names=None):
        super().__init__(properties, hierarchy)
        self.names = names if names is not None else {}
        self._depth_cache = {}

    def __repr__(self) -> str:
        result = super().__repr__()
        return result.replace("SceneParameters", self.__class__.__name__)

    def copy(self):
        return EradiateSceneParameters(
            dict(self.properties),
            {k: list(v) for k, v in self.hierarchy.items()},
            names=dict(self.names),
        )

    def _depth(self, node):
        """
        Recursively compute the maximum depth of a node by walking through all
        its parents. Cache the result in ``_depth_cache`` to amortize cost.
        """
        if node not in self._depth_cache:
            parents = self.hierarchy.get(node, [])
            self._depth_cache[node] = 1 + max(
                (self._depth(p) for p, _ in parents if p is not None), default=-1
            )
        return self._depth_cache[node]

    def set_dirty(self, key: str, expanded=None):
        value, _, node, flags = self.properties[key]

        is_nondifferentiable = flags & mi.ParamFlags.NonDifferentiable
        if is_nondifferentiable and dr.grad_enabled(value):
            mi.Log(
                mi.LogLevel.Warn,
                f"Parameter '{key}' is marked as non-differentiable but has "
                "gradients enabled, unexpected results may occur!",
            )

        self.nodes_to_update.setdefault(node, set()).add(key.rsplit(".", 1)[-1])

        # Climb the tree through all parent nodes using a stack based approach.
        stack = list(self.hierarchy[node])
        if expanded is None:
            expanded = set()
        while stack:
            parent, name = stack.pop()
            if parent is None:
                continue
            self.nodes_to_update.setdefault(parent, set()).add(name)
            if parent in expanded:
                continue
            expanded.add(parent)
            stack.extend(self.hierarchy[parent])

        return self.properties[key]

    def update(self, values=None):
        if values is not None:
            for k, v in values.items():
                if k in self:
                    self[k] = v

        update_candidate_keys = list(self.update_candidates.keys())
        expanded = set()
        for key in update_candidate_keys:
            # Candidate objects might have been modified inplace, we must check
            # the JIT identifiers to see if the object has truly changed.
            if _jit_id_hash(self[key]) == self.update_candidates[key]:
                continue

            self.set_dirty(key, expanded)

        for key in self.keys():
            # explicitely reproduce __get_value__ as it cannot be called.
            value, value_type, node, _ = self.properties[key]
            if value_type is not None:
                value = self.get_property(value, value_type, node)
            dr.schedule(value)

        work_list = sorted(
            self.nodes_to_update.items(),
            key=lambda item: self._depth(item[0]),
            reverse=True,
        )
        out = []
        for node, keys in work_list:
            node.parameters_changed(list(keys))
            out.append((node, keys))

        self.nodes_to_update.clear()
        self.update_candidates.clear()
        dr.eval()

        return out


def traverse(node: mi.Object, name_id_override=None) -> EradiateSceneParameters:
    """
    Like :py:func:`mitsuba.traverse`, but returns an
    :py:class:`EradiateSceneParameters`: a subsequent call to
    :py:meth:`~EradiateSceneParameters.update` notifies every parent of a
    shared object, not just the first one traversal happened to discover.

    Parameter ``node`` (``mi.Object``):
        Node to be traversed.

    Parameter ``name_id_override`` (``bool``, ``str``, ``[str]``, ``None``):
        If set, this argument will be used to select nodes in the scene tree
        whose names will be "pinned" to their ID. Passed values are used as
        regular expressions, with all that it implies regarding ID string
        matching. If this parameter is set to ``True``, a regex that matches
        anything is used.

    Returns → :py:class:`EradiateSceneParameters`:
        The scene parameters obtained from traversing ``node``.
    """

    import re

    if name_id_override is None or name_id_override is False:
        name_id_override = []

    if name_id_override is True:
        name_id_override = [r".*"]

    if type(name_id_override) is not list:
        name_id_override = [name_id_override]

    regexps = [re.compile(k).match for k in name_id_override]

    class EradiateSceneTraversal(mi.TraversalCallback):
        def __init__(
            self,
            node,
            parent=None,
            properties=None,
            hierarchy=None,
            prefixes=None,
            name=None,
            local_name=None,
            flags=+mi.ParamFlags.Differentiable,
            names=None,
        ):
            mi.TraversalCallback.__init__(self)
            self.properties = dict() if properties is None else properties
            self.hierarchy = dict() if hierarchy is None else hierarchy
            self.prefixes = set() if prefixes is None else prefixes
            self.names = dict() if names is None else names

            node_id = node.id()
            if name_id_override and node_id:
                for r in regexps:
                    if r(node_id):
                        name = node_id
                        break

            if name is not None:
                ctr, name_len = 1, len(name)
                while name in self.prefixes:
                    name = "%s_%i" % (name[:name_len], ctr)
                    ctr += 1
                self.prefixes.add(name)

            self.name = name
            self.node = node
            # fully resolved (deduplicated) name.
            self.names[node] = name
            # An object can be reached from more than one parent (shared via
            # a scene-level reference): record every parent edge, keyed by
            # the name under which that specific parent refers to it.
            self.hierarchy.setdefault(node, []).append((parent, local_name))
            self.flags = flags

        def put(self, name, value, flags, cpptype=None):
            if isinstance(value, mi.Object):
                self.put_object(name, value, flags)
            else:
                self.put_value(name, value, flags, cpptype)

        def put_value(self, name, ptr, flags, cpptype):
            name = name if self.name is None else self.name + "." + name

            flags = self.flags | flags
            if (flags & mi.ParamFlags.NonDifferentiable) != 0:
                flags = flags & ~mi.ParamFlags.Discontinuous

            self.properties[name] = (ptr, cpptype, self.node, self.flags | flags)

        def put_object(self, name, obj, flags):
            if obj is None:
                return
            if obj in self.hierarchy:
                # Exists in hierarchy: add parent edge but don't traverse again.
                self.hierarchy[obj].append((self.node, name))
                return
            cb = EradiateSceneTraversal(
                node=obj,
                parent=self.node,
                properties=self.properties,
                hierarchy=self.hierarchy,
                prefixes=self.prefixes,
                name=name if self.name is None else self.name + "." + name,
                local_name=name,
                flags=self.flags | flags,
                names=self.names,
            )
            obj.traverse(cb)

    cb = EradiateSceneTraversal(node)
    node.traverse(cb)

    return EradiateSceneParameters(cb.properties, cb.hierarchy, names=cb.names)
