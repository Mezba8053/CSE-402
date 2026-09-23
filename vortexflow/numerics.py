"""Course numerical methods implemented directly in Python/NumPy."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


def trapezoidal(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("x and y must be one-dimensional arrays of equal length >= 2")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing")
    return float(np.sum(0.5 * np.diff(x) * (y[:-1] + y[1:])))


def simpson_uniform(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != len(y) or len(x) < 3 or (len(x) - 1) % 2:
        raise ValueError("Simpson's rule requires an even number of equal intervals")
    spacing = np.diff(x)
    if not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("Simpson's rule requires uniformly spaced x values")
    return float(
        spacing[0] / 3.0
        * (y[0] + y[-1] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-1:2]))
    )


def lagrange_evaluate(x_nodes: np.ndarray, y_nodes: np.ndarray, x: np.ndarray | float) -> np.ndarray:
    xn = np.asarray(x_nodes, dtype=float)
    yn = np.asarray(y_nodes, dtype=float)
    target = np.asarray(x, dtype=float)
    if len(xn) != len(yn) or len(np.unique(xn)) != len(xn):
        raise ValueError("interpolation nodes must be distinct and match y values")
    result = np.zeros_like(target, dtype=float)
    for i in range(len(xn)):
        basis = np.ones_like(target, dtype=float)
        for j in range(len(xn)):
            if i != j:
                basis *= (target - xn[j]) / (xn[i] - xn[j])
        result += yn[i] * basis
    return result


def newton_coefficients(x_nodes: np.ndarray, y_nodes: np.ndarray) -> np.ndarray:
    x = np.asarray(x_nodes, dtype=float)
    coef = np.asarray(y_nodes, dtype=float).copy()
    if len(x) != len(coef) or len(np.unique(x)) != len(x):
        raise ValueError("interpolation nodes must be distinct and match y values")
    for order in range(1, len(x)):
        coef[order:] = (coef[order:] - coef[order - 1:-1]) / (x[order:] - x[:-order])
    return coef


def newton_evaluate(x_nodes: np.ndarray, coefficients: np.ndarray, x: np.ndarray | float) -> np.ndarray:
    nodes = np.asarray(x_nodes, dtype=float)
    coef = np.asarray(coefficients, dtype=float)
    target = np.asarray(x, dtype=float)
    result = np.full_like(target, coef[-1], dtype=float)
    for i in range(len(coef) - 2, -1, -1):
        result = result * (target - nodes[i]) + coef[i]
    return result


@dataclass(slots=True)
class NaturalCubicSpline:
    x: np.ndarray
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    d: np.ndarray

    def evaluate(self, values: np.ndarray | float) -> np.ndarray:
        values_array = np.asarray(values, dtype=float)
        indices = np.searchsorted(self.x, values_array, side="right") - 1
        indices = np.clip(indices, 0, len(self.a) - 1)
        delta = values_array - self.x[indices]
        return (
            self.a[indices] + self.b[indices] * delta
            + self.c[indices] * delta**2 + self.d[indices] * delta**3
        )


def natural_cubic_spline(x_nodes: np.ndarray, y_nodes: np.ndarray) -> NaturalCubicSpline:
    x = np.asarray(x_nodes, dtype=float)
    y = np.asarray(y_nodes, dtype=float)
    if len(x) != len(y) or len(x) < 3 or np.any(np.diff(x) <= 0):
        raise ValueError("spline nodes must be strictly increasing with at least 3 points")
    n = len(x)
    h = np.diff(x)
    matrix = np.zeros((n, n), dtype=float)
    rhs = np.zeros(n, dtype=float)
    matrix[0, 0] = matrix[-1, -1] = 1.0
    for i in range(1, n - 1):
        matrix[i, i - 1] = h[i - 1]
        matrix[i, i] = 2.0 * (h[i - 1] + h[i])
        matrix[i, i + 1] = h[i]
        rhs[i] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    c_nodes = np.linalg.solve(matrix, rhs)
    a = y[:-1].copy()
    b = np.empty(n - 1)
    d = np.empty(n - 1)
    for i in range(n - 1):
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c_nodes[i] + c_nodes[i + 1]) / 3.0
        d[i] = (c_nodes[i + 1] - c_nodes[i]) / (3.0 * h[i])
    return NaturalCubicSpline(x, a, b, c_nodes[:-1], d)


def qr_least_squares(x: np.ndarray, y: np.ndarray, degree: int) -> tuple[np.ndarray, float]:
    """Polynomial least squares using explicit modified Gram-Schmidt QR."""

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if degree < 0 or degree >= len(x):
        raise ValueError("degree must be non-negative and smaller than sample count")
    design = np.vander(x, N=degree + 1, increasing=True)
    rows, columns = design.shape
    q = np.zeros((rows, columns), dtype=float)
    r = np.zeros((columns, columns), dtype=float)
    work = design.copy()
    for j in range(columns):
        r[j, j] = np.linalg.norm(work[:, j])
        if r[j, j] < 1e-14:
            raise ValueError("rank-deficient least-squares system")
        q[:, j] = work[:, j] / r[j, j]
        for k in range(j + 1, columns):
            r[j, k] = np.dot(q[:, j], work[:, k])
            work[:, k] -= r[j, k] * q[:, j]
    coefficients = np.linalg.solve(r, q.T @ y)
    residual = float(np.linalg.norm(design @ coefficients - y))
    return coefficients, residual


def polynomial_evaluate(coefficients: np.ndarray, x: np.ndarray | float) -> np.ndarray:
    coefficients = np.asarray(coefficients, dtype=float)
    target = np.asarray(x, dtype=float)
    result = np.zeros_like(target, dtype=float)
    for coefficient in coefficients[::-1]:
        result = result * target + coefficient
    return result


def polynomial_derivative(coefficients: np.ndarray, order: int = 1) -> np.ndarray:
    result = np.asarray(coefficients, dtype=float).copy()
    for _ in range(order):
        result = np.array([i * result[i] for i in range(1, len(result))], dtype=float)
        if not len(result):
            return np.array([0.0])
    return result


@dataclass(slots=True)
class OptimizationResult:
    x: float
    value: float
    converged: bool
    iterations: list[dict[str, float]]


def golden_section_maximize(
    function: Callable[[float], float], lower: float, upper: float,
    tolerance: float = 1e-6, max_iterations: int = 200,
) -> OptimizationResult:
    if not lower < upper:
        raise ValueError("lower bound must be below upper bound")
    ratio = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = float(lower), float(upper)
    c = b - ratio * (b - a)
    d = a + ratio * (b - a)
    fc, fd = function(c), function(d)
    history: list[dict[str, float]] = []
    for iteration in range(max_iterations):
        history.append({"iteration": iteration, "lower": a, "upper": b, "c": c, "d": d})
        if abs(b - a) <= tolerance:
            x = 0.5 * (a + b)
            return OptimizationResult(x, function(x), True, history)
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - ratio * (b - a)
            fc = function(c)
        else:
            a, c, fc = c, d, fd
            d = a + ratio * (b - a)
            fd = function(d)
    x = 0.5 * (a + b)
    return OptimizationResult(x, function(x), False, history)


def newton_maximize_polynomial(
    coefficients: np.ndarray, initial: float, lower: float, upper: float,
    tolerance: float = 1e-8, max_iterations: int = 50,
) -> OptimizationResult:
    first = polynomial_derivative(coefficients, 1)
    second = polynomial_derivative(coefficients, 2)
    x = float(initial)
    history: list[dict[str, float]] = []
    for iteration in range(max_iterations):
        gradient = float(polynomial_evaluate(first, x))
        curvature = float(polynomial_evaluate(second, x))
        history.append(
            {"iteration": iteration, "x": x, "gradient": gradient, "curvature": curvature}
        )
        if abs(gradient) <= tolerance and curvature < 0:
            return OptimizationResult(x, float(polynomial_evaluate(coefficients, x)), True, history)
        if abs(curvature) < 1e-14:
            break
        candidate = x - gradient / curvature
        if not lower <= candidate <= upper:
            candidate = min(max(candidate, lower), upper)
            candidate = 0.5 * (x + candidate)
        if abs(candidate - x) <= tolerance:
            x = candidate
            return OptimizationResult(x, float(polynomial_evaluate(coefficients, x)), curvature < 0, history)
        x = candidate
    return OptimizationResult(x, float(polynomial_evaluate(coefficients, x)), False, history)
