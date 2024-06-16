import math
import matplotlib.pyplot as plt
from collections import defaultdict

def plot_points_and_paths(points, paths):
    x = [p[0] for p in points]
    y = [p[1] for p in points]

    _, ax = plt.subplots()

    ax.scatter(x, y, color='blue', label='Points')

    for i, path in enumerate(paths):
        path_x = [p[0] for p in path]
        path_y = [p[1] for p in path]
        path_x.append(path[0][0])
        path_y.append(path[0][1])
        ax.plot(path_x, path_y, label=f'Path {i+1}')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Longest Circular Paths')

    plt.tight_layout()
    plt.show()

def plot_polygon(points, polygon):
    if not points:
        return
    
    x, y = zip(*points)
    
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, color='blue', label='Points')
    
    if polygon:
        px, py = zip(*polygon)
        px += (px[0],)
        py += (py[0],)
        plt.plot(px, py, color='red', label='Convex Hull')
    
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Convex Hull Polygon')
    plt.grid(True)
    plt.show()

def plot_points_and_polygons(points, polygons):
    x = [p[0] for p in points]
    y = [p[1] for p in points]

    _, ax = plt.subplots()

    ax.scatter(x, y, color='blue', label='Points')

    colors = ['red', 'green', 'orange', 'purple', 'cyan']  

    for i, polygon in enumerate(polygons):
        path_x = [p[0] for p in polygon]
        path_y = [p[1] for p in polygon]
        path_x.append(polygon[0][0])  
        path_y.append(polygon[0][1])  
        ax.fill(path_x, path_y, color=colors[i % len(colors)], alpha=0.5, label=f'Polygon {i+1}')  

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Longest Circular Paths and Convex Hull Polygons')

    plt.tight_layout()
    plt.show()

def polar_angle(p0, p1):
    return math.atan2(p1[1] - p0[1], p1[0] - p0[0])

def distance(p0, p1):
    return (p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2

def is_adjacent(p1, p2):
    return distance(p1, p2) == 1 or distance(p1, p2) == 2

def find_bottom_left_point(points):
    return min(points, key=lambda p: (p[1], p[0]))

def sort_points_counterclockwise(points, start_point):
    sorted_points = sorted(points, key=lambda p: polar_angle(start_point, p))
    return sorted_points

def dfs(point, graph, path, visited, start_point, longest_path, global_visited):
    path.append(point)
    visited.add(point)

    for neighbor in graph[point]:
        if neighbor == start_point and len(path) > 2:
            if len(path) > len(longest_path[0]):
                longest_path[0] = path.copy()
        elif neighbor not in visited and neighbor not in global_visited:
            dfs(neighbor, graph, path, visited, start_point, longest_path, global_visited)

    path.pop()
    visited.remove(point)

def build_graph(points):
    graph = defaultdict(list)
    point_set = set(tuple(p) for p in points)

    for point in points:
        point = tuple(point)
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]:
            neighbor = (point[0] + dx, point[1] + dy)
            if neighbor in point_set:
                graph[point].append(neighbor)
    return graph

def find_all_longest_circular_paths(points):
    graph = build_graph(points)
    all_paths = {}
    global_visited = set()

    for start_point in points:
        if start_point not in global_visited:
            visited = set()
            path = []
            longest_path = [[]]
            dfs(start_point, graph, path, visited, start_point, longest_path, global_visited)
            if longest_path[0]:
                sorted_path = longest_path[0]
                all_paths[start_point] = sorted_path
                global_visited.update(longest_path[0])

    return all_paths

def is_valid_segment(p1, p2):
    return (abs(p1[0] - p2[0]) == 1 and abs(p1[1] - p2[1]) == 0) or \
           (abs(p1[0] - p2[0]) == 0 and abs(p1[1] - p2[1]) == 1) or \
           (abs(p1[0] - p2[0]) == 1 and abs(p1[1] - p2[1]) == 1)

def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

def build_hull(points):
    hull = []
    for p in points:
        while len(hull) >= 2 and cross(hull[-2], hull[-1], p) > 0:
            if is_valid_segment(hull[-2], p):
                print("pop")
                hull.pop()
            else:
                break
        hull.append(p)
    return hull

def graham_scan(points):
    return build_hull(points)

POINTS = [
    (4, 11), (5, 10), (5, 11), (5, 12), (6, 10),
    (6, 11), (8, 12), (11, 11), (12, 10), (12, 12),
    (13, 10), (13, 11), (16, 9), (17, 7), (17, 8),
    (17, 10), (18, 8), (18, 9), (19, 10), (21, 9),
    (21, 10), (21, 11), (22, 10), (22, 11),
]

all_paths = find_all_longest_circular_paths(POINTS)
polys = []
for start_point, path in all_paths.items():
    # plot_points_and_paths(POINTS, [path])
    print(path)
    convex_hull_polygon = graham_scan(path)
    polys.append(convex_hull_polygon)
    print(convex_hull_polygon)
    # plot_polygon(POINTS, convex_hull_polygon)

plot_points_and_polygons(POINTS, polys)

