import networkx as nx
import numpy as np
import os
import re
import random
import csv
import subprocess
# TODO S-NODF, W-NODF


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
    return NODF


def null_model(network, source_node, samples=1000):
    """Compute the average and standard deviation of nestedness across an
    ensemble of random replicates within which the interactions of the node
    have been randomised.

    Parameters
    ----------
    network : NetworkX graph

    source_node : node
        The interactions of this node will be randomised.

    samples : int, optional (default = 1000)
        Number of replicates to generate.

    Returns
    -------
    N_mean : float
        Average nestedness of the ensemble.

    N_std : float
        Standard deviation of nestedness of the ensemble.

    References
    ----------
    ..  [1] Bascompte, J., Jordano, P., Melián, C. J. & Olesen, J. M.
        The nested assembly of plant-animal mutualistic networks.
        Proc. Natl Acad. Sci. USA 100, 9383–9387 (2003)
        https://doi.org/10.1073/pnas.1633576100
    """
    N = np.empty(samples)
    same = [node for node in network if network.nodes[node]
            ['bipartite'] == network.nodes[source_node]['bipartite']]
    others = [node for node in network if network.nodes[node]
              ['bipartite'] != network.nodes[source_node]['bipartite']]

    for i in range(samples):
        sample_network = network.copy()
        sample_network.remove_edges_from(network.edges(source_node))
        for other in others:
            p = (network.degree(source_node)/len(others) +
                 network.degree(other)/len(same))/2
            if random.random() < p:
                sample_network.add_edge(source_node, other)
        zero_degree_nodes = [
            node for node in sample_network if sample_network.degree(node) == 0]
        sample_network.remove_nodes_from(zero_degree_nodes)
        N[i] = NODF(sample_network)

    return np.mean(N), np.std(N)


if __name__ == '__main__':
    bin_path = os.path.abspath('bin')
    input_path = os.path.abspath('input')
    output_path = os.path.abspath('output')
    env = os.environ.copy()
    env['PATH'] = '.:' + env['PATH']
    count = 0

    for filename in os.listdir(input_path):
        if (re.match('^.*\.csv$', filename, re.IGNORECASE)):
            print(filename)
            network = nx.Graph()
            lenders = []
            borrowers = []
            loans = []
            with open(os.path.join(input_path, filename), 'r') as edgelist:
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

            sorted_banks = sorted([node for node, data in network.nodes(
                data=True) if data["bipartite"] == 0], key=network.degree, reverse=True)
            sorted_firms = sorted([node for node, data in network.nodes(
                data=True) if data["bipartite"] == 1], key=network.degree, reverse=True)
            A = nx.bipartite.biadjacency_matrix(
                network, row_order=sorted_banks, column_order=sorted_firms).toarray()
            np.savetxt(os.path.join(bin_path, 'matrix' +
                       str(count) + '.in.txt'), A, fmt='%d')
            subprocess.run('contributions -i ' + str(count) +
                           ' --contributions --null-model 0', cwd=bin_path, shell=True, env=env)
            subprocess.run('contributions -i ' + str(count) +
                           ' --nodf', cwd=bin_path, shell=True, env=env)
            row_contributions = []
            column_contributions = []
            c = []
            with open(os.path.join(bin_path, 'matrix' + str(count) + '.contributions.csv'), 'r', newline='') as input_file:
                csv_reader = csv.reader(input_file, delimiter=',')
                next(csv_reader)
                for row in csv_reader:
                    if (row[0] == 'row'):
                        row_contributions.append(row[1])
                    else:
                        column_contributions.append(row[1])
            for id, contribution in zip(sorted_banks, row_contributions):
                c.append({'id': id, 'contribution': contribution, 'type': 'lender'})
            for id, contribution in zip(sorted_firms, column_contributions):
                c.append({'id': id, 'contribution': contribution, 'type': 'borrower'})
            with open(os.path.join(bin_path, 'matrix' + str(count) + '.nodf.txt'), 'r') as input_file:
                nodf = input_file.readline()

            output_folder = os.path.join(
                output_path, os.path.splitext(filename)[0])
            os.makedirs(output_folder, exist_ok=True)
            for output_filename in os.listdir(bin_path):
                if (re.match('^matrix' + str(count) + '.*', output_filename, re.IGNORECASE)):
                    try:
                        os.remove(os.path.join(bin_path, output_filename))
                    except FileNotFoundError:
                        pass

            contributions_filename = os.path.join(
                output_folder, 'contributions.csv')
            with open(contributions_filename, 'w', newline='') as contributions_file:
                csv_writer = csv.DictWriter(contributions_file, fieldnames=[
                                            'id', 'contribution', 'type'])
                csv_writer.writeheader()
                for contribution in c:
                    csv_writer.writerow(contribution)

            nodf_filename = os.path.join(
                output_folder, 'nodf.txt')
            with open(nodf_filename, 'w') as nodf_file:
                nodf_file.write(nodf)
