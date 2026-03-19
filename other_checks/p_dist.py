from simRec import (
    _poisson_draw,
)

lam=1

for _ in range(500):
    print(_poisson_draw(lam))
