def in_range(in_xy,txy):
    return ((in_xy >= min(txy)) & (in_xy <= max(txy)))


def weighted_mean(data, dims, weights):
    weight_sum = weights.sum(dim=dims) # to avoid dividing by zero
    return (data*weights).sum(dim=dims)/weight_sum.where(weight_sum != 0)