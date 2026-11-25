
# SPEC — uniform embedding for Robinson matrices (`robinson-unif-emb`)

This is the 0th version of a specification file, which you will update and improve as the steps progress. Save your next version to `SPEC.md`. 

You are helping the owner transfer one of their papers, `Janssen–Zhang, arXiv-2105.09197v1`, which proposed an algorithm for finding a "uniform embedding for Robinson matrices", into a Python implementation. 

Your first task is to analyze the paper, whose LaTeX source is at `./sources/arXiv-2105.09197v1/ArXiv_ver.tex`. Create documentation in the `./docs` folder and extract all mathematical terms for future reference, either for readers or other language models. In this part, include mathematical definitions and pseudocode that mimics a Python implementation. 

After doing so, do not write any code yet: write the design and interface of the algorithm(s) first. This will be part of your final `SPEC.md` output so I can feed it to the next chat for implementation. 

## Layers 

### Algorithm Core

Here, assume a perfectly suitable input is passed into the function/methods; the pseudocode algorithm is written, and its correctness is verified in the original source, so focus only on the future implementation (inputs, outputs, and the implementation). The transition should be straightforward and focus on the algorithm rather than software engineering; use minimal try-catch-exception in the core. 

### Outer layer

This is the layer where users may use Jupyter Notebooks to code matrices and other vectors (arrays) into the inner layers. It should be simple, input-focused, and friendly to use. 

### Intermediate layer (optional)

The user-friendly requirement and the core input may not coincide. For instance, in another project I encountered, the graph objects need not be complex NetworkX objects but merely edge lists; hence, the intermediate layer in that project translates an `nx` object into a list of edge lists to feed into the algorithm core. 

## Naming

In the paper, the names may not be the most suitable or follow any Python naming convention. Your task, alongside the API design, is also to design appropriate names for methods and functions. This project will be essentially uploaded to GitHub, so in the naming, it should implicitly make suggestions if the function is designed for the project itself but not for general using. 

[The end of v0 spec]