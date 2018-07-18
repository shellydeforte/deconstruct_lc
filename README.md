# deconstruct_lc
Scripts associated with the paper "Deconstructing protein sequence complexity to reveal the building blocks of structure and phase separation"

This project is currently in development.

# Experiments

The paper "Widespread reorganization of metabolic enzymes into reversible assemblies upon nutrient starvation" by Narayanaswamy, Marcotte et al. is the basis for our experiments. In this paper, they started with an undefined group of metabolic yeast proteins (we need to contact them to see if we can get the original list) that consisted of approximately 800 proteins. They observed the formation of puncta after 48 hours of glucose starvation in 180 of these proteins. 

When we apply the low complexity score to these proteins, we can see that the puncta forming proteins are significantly upshifted compared to the yeast proteome. However, there are still a substantial amount that fall below a score of zero. The purpose of these experiments is to try to find out if these lower scoring protein are somehow different. Therefore, starting with the original 180 proteins, we will perform the following experiments.

1. We will check for the formation of puncta before starvation, and we will verify Marcotte's results after starvation.
2. We will check for ThT staining
3. For those proteins that do not stain with ThT, we will check for dissolution under hexandiol treatment

The hypothesis is that high scoring proteins will be more likely to form puncta without stress, less likely to stain with ThT and more likely to dissolve with hexandiol.

Here are the current to do items:
1. Create a spreadsheet with the marcotte proteins and scores
2. Get a list of metabolic proteins for the background
3. Perform chi-square or fisher's exact (probably fisher's exact) on each level of data
4. Re-read the aggregate paper to check for their methods
