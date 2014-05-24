def coerce_bigoh(left,right,op):
    parent = left.parent()
    path = [ ]
    try:
        # Todo: find something better
        path = parent._models.shortest_path(left.__class__,right.__class__)
        if len(path) > 0:
            for i in range(1,len(path)):
                map = parent._models.edge_label(path[i-1],path[i])
                if map == None:
                    left = path[i](parent,left)
                else:
                    left = map(parent,left)
            return op(left,right)
        path = parent._models.shortest_path(right.__class__,left.__class__)
        if len(path) > 0:
            for i in range(1,len(path)):
                map = parent._models.edge_label(path[i-1],path[i])
                if map == None:
                    right = path[i](parent,right)
                else:
                    right = map(parent,right)
            return op(left,right)
    except RuntimeError:
        pass
    raise TypeError("bigoh coercion error: models are %s and %s" % (left._repr_model(), right._repr_model()))
