from parameter import *

'''This script is based on index-0'''

def inverse_mapping(i, Nx):
     
    x = (i) % Nx    
    y = (i) // Nx    
    return (x, y)
 

def neighbor_square(x, y, Nx, Ny, periodic=True):
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

    return sites

print(neighbor_square(2, 2, 5, 4))

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

print(second_nearest(0,0,5,4))