import math

print("md ice 0.001 0.5 0.5 1.0 0.5")
size = 1
vertices = []
for x in range(-size, size+1):
    for y in range(-size, size+1):
        for z in range(-size, size+1):
            fx = float(x)/size
            fy = float(y)/size
            fz = float(z)/size
            print("v", fx/10.0, -fz/10.0, fy/10.0)
            vertices.append((fx, fy, fz))

def normal(a, b, c):
    v = tuple(x-y for x, y in zip(a, b))
    w = tuple(x-y for x, y in zip(c, a))
    return (
        v[1]*w[2]-v[2]*w[1],
        v[2]*w[0]-v[0]*w[2],
        v[1]*w[1]-v[1]*w[0]
    )

def distance(a, b):
    ix, iy, iz = a
    jx, jy, jz = b
    return math.sqrt(pow(ix-jx, 2) + pow(iy-jy, 2) + pow(iz-jz, 2))

def add_triangle(pairs, a, b, c):
    try_add_ordered_triangle(pairs, a, b, c) # or try_add_ordered_triangle(pairs, b, c, a)

def try_add_ordered_triangle(pairs, a, b, c):
    if (a, b) not in pairs and (b, c) not in pairs and (c, a) not in pairs:
        pairs.append((a, b))
        pairs.append((b, c))
        pairs.append((c, a))
        print("f",a+1,b+1,c+1)
        return True
    return False

pairs = []
closest = 100000
for i in range(len(vertices)):
    for j in range(len(vertices)):
        for k in range(len(vertices)):
            distances = [
                distance(vertices[i], vertices[j]),
                distance(vertices[j], vertices[k]),
                distance(vertices[k], vertices[i])
            ]
            distances.sort()
            if distances[0] == 1.0 and distances[1] == 1.0 and distances[2] == math.sqrt(2):
                n = normal(vertices[i], vertices[j], vertices[k])
                # print(n)
                if (vertices[i][1] == 1.0 and vertices[j][1] == 1.0 and vertices[k][1] == 1.0
                    and n[0] == 0.0 and n[1] == -1.0 and n[2] == 0.0):
                    add_triangle(pairs, i, j, k)
                elif (vertices[i][1] == -1.0 and vertices[j][1] == -1.0 and vertices[k][1] == -1.0
                    and n[0] == 0.0 and n[1] == 1.0 and n[2] == 0.0):
                    add_triangle(pairs, i, j, k)
                elif (vertices[i][0] == 1.0 and vertices[j][0] == 1.0 and vertices[k][0] == 1.0
                    and n[0] == -1.0 and n[1] == 0.0 and n[2] == 0.0):
                    add_triangle(pairs, i, j, k)
                elif (vertices[i][0] == -1.0 and vertices[j][0] == -1.0 and vertices[k][0] == -1.0
                    and n[0] == 1.0 and n[1] == 0.0 and n[2] == 0.0):
                    add_triangle(pairs, i, j, k)
                # elif (vertices[i][2] == 1.0 and vertices[j][2] == 1.0 and vertices[k][2] == 1.0
                #     and n[0] == 0.0 and n[1] == 0.0 and n[2] == -1.0):
                #     add_triangle(pairs, i, j, k)
                elif (vertices[i][2] == -1.0 and vertices[j][2] == -1.0 and vertices[k][2] == -1.0
                    and n[0] == 0.0 and n[1] == 0.0 and n[2] == 1.0):
                    add_triangle(pairs, i, j, k)
                # break