# 不需要特殊配置，直接导入 Sage 库
from sage.all import *

# 使用 Sage 功能
x = var('x')
expr = sin(x**2) + cos(x)
print(f"表达式: {expr}")
print(f"导数: {expr.diff(x)}")

# 符号积分
integral_result = integral(exp(-x**2), (x, -oo, oo))
print(f"高斯积分: {integral_result}")
print(f"数值近似: {integral_result.n()}")
