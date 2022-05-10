# Fetch and import APX wrapper class
import apx
from importlib import reload

reload(apx)

from typing import Iterable
from apx import DataFile, LinearProgram, np


###############################################################################
# file reader function

def file_reader(filename: str) -> list:
    """Reads a given file.

    Args:
        filename (str): path and name of the file.

    Returns:
        list(list(str)): lines of the file, where each line is separated in 
            sections by the whitespace.
    """
    with open(filename, 'r', encoding='utf8') as datafile:
        data = datafile.read().splitlines()
    
    return [line.split(" ") for line in data]


###############################################################################
# Implementing the Min Buckets algorithm (trivial algorithm)

def find_degrees(edges: Iterable) -> dict:
    """Finds the degree of the vertices of a graph.

    Args:
        edges (Iterable): pairs of vertices of the graph that correspond to an
            edge.

    Returns:
        vertices (dict): vertices with associated degree.
    """
    vertices = {}
    for (u, v) in edges:
        if u in vertices:
            vertices[u] += 1
        else:
            vertices[u] = 1
        if v in vertices:
            vertices[u] += 1
        else:
            vertices[v] = 1
    return vertices


def find_wedges(edges: Iterable, vert_deg: dict) -> dict:
    """For every pair of vertices that correspond to an edge in the graph, 
    associate lowest degree vertice to all the higher degree vertices that 
    have an edge with it.

    Args:
        edges (Iterable): pairs of vertices of the edges of the graph.
        vert_deg (dict): vertices with their associated degrees.

    Returns:
        wedges (dict): lowest degree vertices associated to all the highest
            degree vertices with which it forms an edge in the graph.
    """

    wedges = {}
    for (u, v) in edges:
        least, highest = (u,v) if vert_deg[u] < vert_deg[v] else (v, u)
        if least in wedges:
            wedges[least].add(highest)
        else:
            wedges[least] = {highest}
    return wedges


def find_triangles(edges: Iterable, wedges: dict) -> set:
    """Finds all the triangles in a graph.

    Args:
        edges (Iterable): pairs of vertices that form edges in the graph.
        wedges (dict): lowest degree vertices associated to all the highest
            degree vertices with which it forms an edge in the graph.

    Returns:
        triangle (set): triples of vertices that form all triangles (all 
            combinations of two of the three vertices correspond to an edge of
            the graph) of the graph.
    """

    triangle = set()
    for edge in edges:
        for vert in wedges:
            if all([vertice in wedges[vert] for vertice in edge]):
                triangle.add(frozenset((vert, *edge)))
    return triangle


def min_buckets(edges: Iterable) -> set:
    """Applies the Minimum Buckets Algorithm (trivial algorithm) to find all
    instances of triangles in a graph.

    Args:
        edges (Iterable): pairs of vertices that form edges in the graph.

    Returns:
        triangles (set): triples of vertices that form all triangles (all 
            combinations of two of the three vertices correspond to an edge of
            the graph) of the graph.
    """
    vert_deg = find_degrees(edges)
    wedges = find_wedges(edges, vert_deg)
    triangles = find_triangles(edges, wedges)
    return triangles


###############################################################################
# solving Vertex Cover using Linear Programming

def solve_vertex_cover(filename: str) -> tuple:
    """Original implementation of the solution of the Vertex Cover problem, 
    using linear programming."""
    
    graph = DataFile(filename)
    vertex_cover_lp = LinearProgram('min')
    objective = {}
    for (u,v) in graph:
        vertex_cover_lp.add_constraint({u: 1, v: 1}, 1)
        objective[u] = 1.0
        objective[v] = 1.0
    
    vertex_cover_lp.set_objective(objective)
    value, solution = vertex_cover_lp.solve()
    return value, solution


def solve_vertex_cover_triangles(filename: str) -> tuple:
    """Implementation of the solution of the Vertex Cover problem, 
    using linear programming, with the additional constraint that two vertices
    of the three vertices that form a triangle in the graph must be selected, 
    for each triangle."""
    
    edges = file_reader("data/" + filename)
    graph = DataFile(filename)
    vertex_cover_lp = LinearProgram('min')
    objective = {}
    for (u,v) in graph:
        vertex_cover_lp.add_constraint({u: 1, v: 1}, 1)
        objective[u] = 1.0
        objective[v] = 1.0

    for triangle in min_buckets(edges):
        vertex_cover_lp.add_constraint({edge: -1 for edge in triangle}, -2)
    
    vertex_cover_lp.set_objective(objective)
    value, solution = vertex_cover_lp.solve()
    return value, solution


def rounding(solution) -> tuple:
    """Performs the rounding of a solution to the vertex cover problem."""
    
    rounded_value, rounded_solution = 0, {}
    for x in solution:
        r = int(np.round(solution[x] + 1e-10)) 
        # Add small constant to deal with numerical issues for numbers close to 1/2
        rounded_solution[x] = r
        rounded_value += r
    return rounded_value, rounded_solution


###############################################################################

if __name__ == '__main__':
    for filename in DataFile.graph_files:
        print(f'Graph: {filename}')
        val, sol = solve_vertex_cover(filename)
        round_val, round_sol = rounding(sol)

        t_val, t_sol = solve_vertex_cover_triangles(filename)
        t_round_val, t_round_sol = rounding(t_sol)
        
        print(f"")
        print(f'LP optimal value: original -> {val}; triangle -> {t_val}')
        print('original')
        print(f'{sol}')
        print('triangle')
        print(f'{t_sol}')
        print()
        print(f'Rounded LP value: original -> {round_val}; triangle ->' +
                f'{t_round_val}')
        print('original')
        print(round_sol)
        print('triangle')
        print(t_round_sol)
        print()
        print()
        print("-"*20)