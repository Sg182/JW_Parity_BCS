from parameter import *

'''This script is based on index-0'''

def inverse_mapping(i, Nx):
     
    x = (i) % Nx    
    y = (i) // Nx    
    return (x, y)
 

'''def neighbor_square(x, y, Nx, Ny, periodic=True):
    sites = []

    # Map (x,y) → index i
    def map_square(x, y, Nx,Ny):
        assert y <Ny, f"y={y} out of bounds"
        return y * Nx + x

    if periodic:
        # Right neighbor
        nx = (x + 1) % Nx
        sites.append(map_square(nx, y, Nx,Ny))

        # Left neighbor
        nx = (x - 1) % Nx
        sites.append(map_square(nx, y, Nx,Ny))

        # Up neighbor
        ny = (y + 1) % Ny
        sites.append(map_square(x, ny, Nx,Ny))

        # Down neighbor
        ny = (y - 1) % Ny
        sites.append(map_square(x, ny, Nx,Ny))

    else:
        # Right neighbor
        if x + 1 < Nx:
            sites.append(map_square(x + 1, y, Nx,Ny))

        # Left neighbor
        if x - 1 >= 0:
            sites.append(map_square(x - 1, y, Nx,Ny))

        # Up neighbor
        if y + 1 < Ny:
            sites.append(map_square(x, y + 1, Nx,Ny))

        # Down neighbor
        if y - 1 >= 0:
            sites.append(map_square(x, y - 1, Nx,Ny))

    return sites'''

#print(neighbor_square(2, 2, 5, 4))

def neighbors_square_1st(x, y, Nx, Ny, periodic=True):
    # ---- Validate INPUT coordinates ----
    '''For 1D set Nx or Ny to 1'''
    assert 0 <= x < Nx, f"input x={x} out of bounds (Nx={Nx})"
    assert 0 <= y < Ny, f"input y={y} out of bounds (Ny={Ny})"

    nn = []

    def map_square(x, y, Nx, Ny):
        return y * Nx + x

    def wrap_x(u): return (u + Nx) % Nx if periodic else u
    def wrap_y(v): return (v + Ny) % Ny if periodic else v

    # Horizontal neighbors
    if Nx > 1:
        for nx in (x+1, x-1):
            if periodic or (0 <= nx < Nx):
                nn.append(map_square(wrap_x(nx), y, Nx, Ny))

    # Vertical neighbors
    if Ny > 1:
        for ny in (y-1, y+1):
            if periodic or (0 <= ny < Ny):
                nn.append(map_square(x, wrap_y(ny), Nx, Ny))

    return nn

print(neighbors_square_1st(0,1,5,2))



def second_nearest(x, y, Nx, Ny, periodic=False):

    '''returns second-nearest points for 
    a bipartitie lattice.'''
    diag = []
    def map_square(x, y, Nx, Ny):
        assert 0 <= x < Nx, f"x={x} out of bounds"
        assert 0 <= y < Ny, f"y={y} out of bounds"
        return y * Nx + x

    def wrap_x(u): return (u + Nx) % Nx if periodic else u
    def wrap_y(v): return (v + Ny) % Ny if periodic else v

    candidates = [(x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1)]
    for nx, ny in candidates:
        if periodic:
            diag.append(map_square(wrap_x(nx), wrap_y(ny), Nx, Ny))
        else:
            if 0 <= nx < Nx and 0 <= ny < Ny:
                diag.append(map_square(nx, ny, Nx, Ny))
    return diag


def neighbors_2nd_square(x, y, Nx, Ny, periodic=False):
    """
    2D: diagonals (x±1, y±1).
    1D (Ny==1): x±2 along x.
    1D (Nx==1): y±2 along y.
    FOR 1D, set Nx or Ny to 1
    """
    assert Nx >= 1 and Ny >= 1, f"Invalid lattice size Nx={Nx}, Ny={Ny}"
    assert 0 <= x < Nx and 0 <= y < Ny, f"(x,y)=({x},{y}) out of bounds"
    def map_square(x, y, Nx, Ny):
        assert 0 <= x < Nx, f"x={x} out of bounds [0,{Nx-1}]"
        assert 0 <= y < Ny, f"y={y} out of bounds [0,{Ny-1}]"
        return y * Nx + x

    out = []

    def wrap_x(u): return (u + Nx) % Nx if periodic else u
    def wrap_y(v): return (v + Ny) % Ny if periodic else v

    if Nx > 1 and Ny > 1:
        # 2D: diagonals
        candidates = [(x+1, y+1), (x+1, y-1), (x-1, y+1), (x-1, y-1)]
        for nx, ny in candidates:
            if periodic:
                out.append(map_square(wrap_x(nx), wrap_y(ny), Nx, Ny))
            else:
                if 0 <= nx < Nx and 0 <= ny < Ny:
                    out.append(map_square(nx, ny, Nx, Ny))

    elif Ny == 1 and Nx > 1:
        # 1D along x: ±2
        for nx in (x + 2, x - 2):
            if periodic or (0 <= nx < Nx):
                out.append(map_square(wrap_x(nx), 0, Nx, 1))

    elif Nx == 1 and Ny > 1:
        # 1D along y: ±2
        for ny in (y + 2, y - 2):
            if periodic or (0 <= ny < Ny):
                out.append(map_square(0, wrap_y(ny), 1, Ny))

    # Nx==Ny==1: single site, no neighbors
    return out


