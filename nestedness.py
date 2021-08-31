import networkx as nx
import numpy as np
import os
import re
import subprocess

input_path = os.path.abspath('input')
bin_path = os.path.abspath('bin')
output_path = os.path.abspath('output')
env = os.environ.copy()
env['PATH'] = '.:' + env['PATH']

for filename in os.listdir(input_path):
    if (re.match('^.*\.csv$', filename, re.IGNORECASE)):
        with open(os.path.join(input_path, filename), 'rb') as edgelist:
            edgelist.readline()  # Skips first line
            network = nx.bipartite.read_edgelist(
                edgelist, delimiter=',', data=[('loan value', float)])

        sorted_banks = sorted([node for node, data in network.nodes(
            data=True) if data["bipartite"] == 0], key=network.degree, reverse=True)
        sorted_firms = sorted([node for node, data in network.nodes(
            data=True) if data["bipartite"] == 1], key=network.degree, reverse=True)
        A = nx.bipartite.biadjacency_matrix(
            network, row_order=sorted_banks, column_order=sorted_firms).toarray()
        np.savetxt(os.path.join(bin_path, 'matrix0ord.txt'), A, fmt='%d')
        np.savetxt(os.path.join(bin_path, 'general0.txt'), A.shape, fmt='%d')

        # Construction of null model
        # TODO use 'matrix0rand.txt' if available
        subprocess.run('simulated_annealing 0 1',
                       cwd=bin_path, shell=True, env=env)

        # Analytical measures
        subprocess.run('NODF_analytic 0', cwd=bin_path, shell=True, env=env)
        subprocess.run('Rscript spectral_radius_analytic.R 0',
                       cwd=bin_path, shell=True, env=env)

        # Sampled measures

        # Low sample size to test computational time required, recommended
        # values are a minimum of 1000 and ideally 10000 for the best balance
        # between accuracy and speed
        Nsample = '100'
        Nreps = '500'
        subprocess.run('Rscript temperature.R 0 ' + Nsample,
                       cwd=bin_path, shell=True, env=env)
        subprocess.run('NODF 0 ' + Nsample, cwd=bin_path, shell=True, env=env)
        subprocess.run('Rscript discrepancy.R 0 ' + Nsample,
                       cwd=bin_path, shell=True, env=env)
        subprocess.run('NIR 0 ' + Nsample, cwd=bin_path, shell=True, env=env)
        subprocess.run('Rscript spectral_radius.R 0 ' + Nsample,
                       cwd=bin_path, shell=True, env=env)
        # Nestedness based on Manhattan distance requires bipartite version 0.90
        subprocess.run('Rscript NMD.R 0 ' + Nsample + ' ' +
                       Nreps, cwd=bin_path, shell=True, env=env)

        output_folder = os.path.join(output_path, os.path.splitext(filename)[0])
        os.makedirs(output_folder, exist_ok=True)
        for output_filename in os.listdir(bin_path):
            if (re.match('^.*\.txt$', output_filename, re.IGNORECASE)):
                os.replace(os.path.join(bin_path, output_filename),
                           os.path.join(output_folder, output_filename))
