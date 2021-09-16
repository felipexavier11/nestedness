import networkx as nx
import numpy as np
import os
import re
import random
import csv
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
    input_path = os.path.abspath('input')
    output_path = os.path.abspath('output')

    for filename in os.listdir(input_path):
        if (re.match('^.*\.csv$', filename, re.IGNORECASE)):
            with open(os.path.join(input_path, filename), 'rb') as edgelist:
                edgelist.readline()  # Skips first line
                # FIXME if firms and banks can have the same ID this will break
                network = nx.bipartite.read_edgelist(
                    edgelist, delimiter=',', data=[('loan value', float)])
            try:
                N = NODF(network)
            except ValueError:
                continue
            print(filename)
            c = []
            for node, data in network.nodes(data=True):
                N_mean, N_std = null_model(network, node, samples=1000)
                node_type = 'bank' if data['bipartite'] == 0 else 'firm'
                c.append({'id': node, 'contribution': (
                    N-N_mean)/N_std, 'type': node_type})

            output_folder = os.path.join(output_path, os.path.splitext(filename)[0])
            os.makedirs(output_folder, exist_ok=True)
            output_file_path = os.path.join(output_folder, 'contributions.csv')
            with open(output_file_path, 'w', newline='') as output_file:
                csv_writer = csv.DictWriter(output_file, fieldnames=[
                                            'id', 'contribution', 'type'])
                csv_writer.writeheader()
                for contribution in c:
                    csv_writer.writerow(contribution)
