"""Indigo independent numerical verification primitives; stdlib only."""
import math

def erdos_step_bound(values, tolerance=1e-8):
    """Feasibility and maximal zero-extended cross-correlation for equal-width steps."""
    h = [float(x) for x in values]
    n = len(h)
    if not n or n % 2 or any(not math.isfinite(x) or x < 0 or x > 1 for x in h):
        raise ValueError("Invalid even-length step function")
    if abs(math.fsum(h) - n / 2) > tolerance:
        raise ValueError("Integral constraint failed")
    best = 0.0
    for shift in range(-(n - 1), n):
        score = math.fsum(h[i] * (1 - h[i + shift]) for i in range(n) if 0 <= i + shift < n)
        best = max(best, score)
    return 2 * best / n

def tammes_min_distance(vectors):
    """Normalize 50 finite nonzero 3-vectors and measure minimum chord distance."""
    if len(vectors) != 50:
        raise ValueError("Expected exactly 50 vectors")
    points = []
    for v in vectors:
        if len(v) != 3 or any(not math.isfinite(float(x)) for x in v):
            raise ValueError("Invalid vector")
        norm = math.sqrt(math.fsum(float(x) ** 2 for x in v))
        if norm == 0:
            raise ValueError("Zero vector")
        points.append(tuple(float(x) / norm for x in v))
    return min(math.dist(a, b) for i, a in enumerate(points) for b in points[i + 1:])

if __name__ == "__main__":
    assert abs(erdos_step_bound([.5, .5]) - .5) < 1e-12
    assert abs(erdos_step_bound([1, 0]) - 1) < 1e-12
    try:
        erdos_step_bound([.3, .3])
        raise AssertionError("Invalid integral accepted")
    except ValueError:
        pass
    try:
        tammes_min_distance([[1, 0, 0]] * 49 + [[0, 0, 0]])
        raise AssertionError("Zero vector accepted")
    except ValueError:
        pass
    print("PASS: four smoke checks")
