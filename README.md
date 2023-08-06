# Binary Comprehensive Learning Particle Swarm Optimization

This **MATLAB** code presents the binary comprehensive learning particle swarm optimization (BCLPSO) method, a powerful meta-heuristic algorithm for the optimal design of nonlinear steel structures using standard member sizes. Meta-heuristic algorithms are a class of advanced optimization techniques that employ trial-and-error processes to search for optimal solutions within a population. The term **meta** signifies a higher level beyond traditional approaches.

The BCLPSO algorithm belongs to the population-based (or trajectory-based) meta-heuristic techniques, characterized by two main phases: exploitation and exploration. The exploitation phase focuses on searching for local solutions by leveraging information from good local solutions, while the exploration phase constructs the search space to discover global optima. The key to effective meta-heuristic methods lies in striking a balance between exploration and exploitation, ensuring the avoidance of premature convergence to local optima and increasing the chances of finding accurate optima.

Numerous meta-heuristic algorithms have been developed, each with its own strengths in terms of exploitation and exploration. Particle swarm optimization (PSO), a swarm-intelligence approach inspired by the collective behavior of bird flocks, is one such algorithm. PSO constructs a population of particles whose positions are iteratively updated based on velocity functions influenced by the global best particle. However, standard PSO often falls prey to premature convergence to local optima due to insufficient social update components. To address this limitation, various enhancements have been proposed to boost the global search capability and overcome local optimal pitfalls.

Our recent work introduces an exceptional variant of PSO, known as comprehensive learning particle swarm optimization (CLPSO), which has demonstrated remarkable performance in various engineering applications. CLPSO incorporates a learning technique that facilitates cross positions between sets of the best swarm particles in each dimensional space, effectively overcoming locally optimal searches and preventing premature termination of non-optimal but feasible solutions. The proposed scheme utilizes a learning probability function to foster cooperative responses among swarm populations, enhancing the algorithm's ability to explore the search space effectively.

By harnessing the power of the BCLPSO algorithm and integrating concepts from the CLPSO variant, this work aims to provide an advanced computational tool for the optimal design of nonlinear steel structures, enabling engineers to efficiently navigate the design space and obtain accurate solutions.

## This code offers several benefits:
- The design complies with the AISC-LRFD standard specifications.
- The sizes and lay-outs of cross-brace members appended to the steel frames are simultaneously optimized.
- The latter converts design variables into binary bit-strings.

## Contributing:
This project welcomes contributions and suggestions from interested individuals. For more detailed information and to delve into the depths of this research, kindly refer to the publication reference provided below:

[1] Su, R.; Tangaramvong, S.; Van, T.H.; Chaiwongnoi, A.; Song, C.  Binary Comprehensive Learning Particle Swam Optimization Approach for Optimal Design of Nonlinear Steel Structures with Standard Sizes. Buildings 2023, 13(8), 1988.
(https://doi.org/10.3390/buildings13081988)
