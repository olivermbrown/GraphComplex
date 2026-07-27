import cells as cls
import configs

N = 10
D11 = cls.TriangleCell(N)
D12 = cls.SquareCell(N)
D22 = cls.TriangleCell(N)
D11.indices = (1, 1)
D12.indices = (1, 2)
D22.indices = (2, 2)
C = configs.ConfigurationSpace([D11, D12, D22])
C.diagonal_boundary_condition = 'neumann'
C.exterior_boundary_condition = 'neumann'

try:
    C.load_states('does/not/exist_', override_directory_path=True)
except Exception as e:
    print(type(e).__name__, e)
