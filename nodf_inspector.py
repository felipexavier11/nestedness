import csv
import os
import re
import subprocess
import networkx as nx
import numpy as np

def NODF(network):
    """Compute the NODF of the network.

    Parameters
    ----------
    network : NetworkX graph

    Returns
    -------
    NODF : float

    Raises
    ------
    ValueError
        If network has zero degree node.

    References
    ----------
    ..  [1] Almeida-Neto, M., Guimarães, P., Guimarães, P.R., Jr, Loyola, R.D.
        and Ulrich, W.
        A consistent metric for nestedness analysis in ecological systems:
        reconciling concept and measurement.
        Oikos, 117: 1227-1239 (2008).
        https://doi.org/10.1111/j.0030-1299.2008.16644.x
    """
    for node in network:
        if network.degree(node) == 0:
            raise ValueError('network has zero degree node')

    row = [node for node, data in network.nodes(
        data=True) if data["bipartite"] == 0]
    col = [node for node, data in network.nodes(
        data=True) if data["bipartite"] == 1]
    n = len(row)
    m = len(col)

    N_row = 0
    for i in range(n-1):
        for j in range(i+1, n):
            # Decreasing fill
            if network.degree(row[i]) != network.degree(row[j]):
                # Paired overlap
                overlap = set(network.neighbors(row[i])) & set(
                    network.neighbors(row[j]))
                k_min = min(network.degree(row[i]), network.degree(row[j]))
                N_row += len(overlap)/k_min

    N_col = 0
    for i in range(m-1):
        for j in range(i+1, m):
            if network.degree(col[i]) != network.degree(col[j]):
                overlap = set(network.neighbors(col[i])) & set(
                    network.neighbors(col[j]))
                k_min = min(network.degree(col[i]), network.degree(col[j]))
                N_col += len(overlap)/k_min

    K = (n*(n-1))/2 + (m*(m-1))/2
    NODF = 100.0*(N_row + N_col)/K
    return N_row, N_col, NODF
    
if __name__ == '__main__':
    bin_path = os.path.abspath('bin')
    input_path = os.path.abspath('input')
    output_path = os.path.abspath('output')
    env = os.environ.copy()
    env['PATH'] = '.:' + env['PATH']
    count = 0

    filename = input('Type the name of the file in the input folder or "exit" (e.g. sample1.csv): ')

    while filename != 'exit':
        with open(os.path.join(input_path, filename), 'r') as edgelist:
            network = nx.Graph()
            lenders = []
            borrowers = []
            loans = []
            edgelist.readline()  # Skips first line
            for edge in edgelist:
                lender, borrower, loan_value = edge.strip().split(',')
                lender = lender + '_lender'
                borrower = borrower + '_borrower'
                lenders.append(lender)
                borrowers.append(borrower)
                loans.append((lender, borrower))
            network.add_nodes_from(lenders, bipartite=0)
            network.add_nodes_from(borrowers, bipartite=1)
            network.add_edges_from(loans)
            lenders_contribution, borrowers_contribution, NODF_value = NODF(network)
            print(f"Lenders' contribution: {lenders_contribution}")
            print(f"Borrowers' contribution: {borrowers_contribution}")
            print(f"NODF for target network: {NODF_value}")
        filename = input('Type the name of the file in the input folder or "exit" (e.g. sample1.csv): ')
