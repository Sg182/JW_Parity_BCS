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
