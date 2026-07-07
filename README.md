# Inference of Multiscale Gaussian Graphical Models
Edmond Sanou, Christophe Ambroise, Geneviève Robin
2023-06-28

### Citation

Edmond Sanou, Christophe Ambroise and Geneviève Robin (June 2023). Inference of Multiscale Gaussian Graphical Models. Computo.
<https://doi.org/10.57750/1f4p-7955>

### Badges

[![build and
publish](https://github.com/computorg/published-202306-sanou-multiscale_glasso/actions/workflows/build.yml/badge.svg)](https://github.com/computorg/published-202306-sanou-multiscale_glasso/actions/workflows/build.yml)
[![reviews](https://img.shields.io/badge/review-report-blue)](https://github.com/computorg/published-202306-sanou-multiscale_glasso/issues?q=is%3Aopen+is%3Aissue+label%3Areview)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/computorg/published-202306-sanou-multiscale_glasso)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/computorg/published-202306-sanou-multiscale_glasso)
[![DOI:10.57750/1f4p-7955](https://img.shields.io/badge/DOI-10.57750%2F1f4p--7955-034E79.svg)](https://doi.org/10.57750/1f4p-7955)
[![Creative Commons
License](https://i.creativecommons.org/l/by/4.0/80x15.png)](http://creativecommons.org/licenses/by/4.0/)

### Authors’ affiliations

- [Edmond Sanou](https://desanou.github.io/) (Université Paris-Saclay, CNRS, Univ Evry, Laboratoire de Mathématiques et Modélisation d’Evry)
- [Christophe Ambroise](https://cambroise.github.io/) (Université Paris-Saclay, CNRS, Univ Evry, Laboratoire de Mathématiques et Modélisation d’Evry)
- [Geneviève Robin](https://genevieverobin.wordpress.com/) (Université Paris-Saclay, CNRS, Univ Evry, Laboratoire de Mathématiques et Modélisation d’Evry)

### Abstract

Gaussian Graphical Models (GGMs) are widely used in high-dimensional
data analysis to synthesize the interaction between variables. In many
applications, such as genomics or image analysis, graphical models rely
on sparsity and clustering to reduce dimensionality and improve
performances. This paper explores a slightly different paradigm where
clustering is not knowledge-driven but performed simultaneously with the
graph inference task. We introduce a novel Multiscale Graphical Lasso
(MGLasso) to improve networks interpretability by proposing graphs at
different granularity levels. The method estimates clusters through a
convex clustering approach — a relaxation of $k$-means, and hierarchical
clustering. The conditional independence graph is simultaneously
inferred through a neighborhood selection scheme for undirected
graphical models. MGLasso extends and generalizes the sparse group fused
lasso problem to undirected graphical models. We use continuation with
Nesterov smoothing in a shrinkage-thresholding algorithm (CONESTA) to
propose a regularization path of solutions along the group fused Lasso
penalty, while the Lasso penalty is kept constant. Extensive experiments
on synthetic data compare the performances of our model to
state-of-the-art clustering methods and network inference models.
Applications to gut microbiome data and poplar’s methylation mixed with
transcriptomic data are presented.
