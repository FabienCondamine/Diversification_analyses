# Diversification_analyses
R codes to tease apart mountain uplift, climate change and biotic drivers of species diversification

Fabien L. Condamine*, Alexandre Antonelli, Laura P. Lagomarsino, Carina Hoorn and Lee Hsiang Liow

Author for correspondence (*): fabien.condamine@gmail.com

Manual to use the R script provided in Dryad
Oct. 15, 2015 – Fabien L. Condamine

The following R script aims at testing the prevalence of the Red Queen and the Court Jester using time-calibrated phylogenies only. The script basically takes a tree (or multiple trees), and performs several birth-death models as explained in the chapter, and detailed below. 

The Court Jester analyses are represented by the PaleoEnv model (Condamine et al., 2013) in which one can set an environmental variable that itself varies through time (e.g. temperature, sea-level). The model will thus estimate whether the speciation/extinction varied according to the past environment, and to which extent.

The Red Queen analyses are made with the diversity-dependent (DDD) model that estimates whether a clade has reached its ecological limits, i.e. its carrying capacity (Etienne et al., 2012). Also, one can easily change the PaleoEnv model into a multiple clade diversity-dependent model if he has past diversity curve of an interacting clade with the focal clade. This will estimate the correlation between the focal clade’s speciation/extinction and the diversity dynamics of the extra clade.

In the book chapter, Condamine et al. (2017) used this general approach, although a simplified one, on the Andean clade of hummingbirds. They tested whether the Andean orogeny, or change in global past temperatures, or the role of intra-clade ecological interactions best explain the diversification of the clade using that framework. 

The framework is highly flexible, as explained in the chapter, and may include other models to test for the prevalence of the Red Queen and Court Jester. For instance, the time-dependent models of Morlon et al. (2011) and Stadler (2011) are included, but were not tested on the Andean clade of hummingbirds. The time-dependent models, PaleoEnv models, and DDD models are comparable based on their corrected Akaike Information Criterion (AICc, see the chapter for an explanation). Note that the TreePar model is not comparable to the others.

Please, cite the book chapter and appropriate references (see below) when using this script and approach. Also, don’t hesitate to contact me if you have any trouble or question.

Etienne R.S., Haegeman B., Stadler T., Aze T., Pearson P.N., Purvis A., Phillimore A.B. 2012. Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record. Proc. Roy. Soc. Lond. B, 279, 1300–1309.
Morlon H., Parsons T.L., Plotkin, J. 2011. Reconciling molecular phylogenies with the fossil record. Proc. Natl. Acad. Sci. USA, 108, 16327-16332.
Condamine F.L., Rolland J., Morlon H. (2013) Macroevolutionary perspectives to environmental change. Ecol. Lett., 16, 72-85.
Condamine F.L., Antonelli A., Lagomarsino L.P., Hoorn C., Liow L.H. 2017. Teasing apart mountain uplift, climate change and biotic drivers of species diversification. In: Mountains, Climate, and Biodiversity (eds. Hoorn C., Antonelli A.). Wiley Blackwell, pp. xx-xx. 
Stadler T. 2011. Mammalian phylogeny reveals recent diversification rate shifts. Proc. Natl. Acad. Sci. USA, 108, 6187-6192.
