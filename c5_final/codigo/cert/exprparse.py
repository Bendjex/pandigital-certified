"""build an exact.X value from a formula string written by motor B (fully parenthesised)"""
import ast, re
from fractions import Fraction as Fr
import exact as ex


def to_x(expr):
    py = re.sub(r"(\d*\.\d+|\d+)", lambda m: f"N('{m.group(1)}')", expr).replace("^", "**")
    def ev(t):
        if isinstance(t, ast.Call):
            return ex.X.rat(Fr(t.args[0].value))
        if isinstance(t, ast.UnaryOp) and isinstance(t.op, ast.USub):
            return ex.neg(ev(t.operand))
        if isinstance(t, ast.BinOp):
            x, y = ev(t.left), ev(t.right)
            if isinstance(t.op, ast.Add): return ex.add(x, y)
            if isinstance(t.op, ast.Sub): return ex.sub(x, y)
            if isinstance(t.op, ast.Mult): return ex.mul(x, y)
            if isinstance(t.op, ast.Div): return ex.div(x, y)
            if isinstance(t.op, ast.Pow): return ex.power(x, y)
        raise ValueError(expr)
    return ev(ast.parse(py, mode="eval").body)
