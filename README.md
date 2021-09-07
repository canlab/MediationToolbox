# MediationToolbox

The multilevel mediation and moderation (M3) toolbox is a Matlab toolbox designed for mediation analysis. It permits tests of both single-level and multi-level mediation. It was developed in part for use with neuroimaging data, but can be used for any type of data, including behavior/performance, survey data, physiological data, and more.

**Features**

The M3 toolbox provides several features tailored for mediation analysis. These are available when analyzing any kind of data, and are also available when performing voxel-wise analyses of neuroimaging data:
- Single or multi-level mediation [mediation.m]
- Bias corrected, accelerated bootstrapping for inference
- Precision-weighted estimates for mixed effects (multi-level)
- Logistic regression for categorical outcomes (Y)
- Ability to add multiple mediators and covariates
- Second-level moderators for moderated mediation
- Multi-path (3-path) mediation
- Permutation testing for inference
- Latent hemodyamic response basis set and time-shift for fMRI mediators
- Autocorrelation estimation (AR(p)) for time series mediators
- Multivariate mediation to identify a pattern across dense, high-dimensional mediators [PDM toolbox](https://github.com/canlab/MediationToolbox/tree/master/PDM_toolbox)

**Main functions to run**

The function mediation.m is the main function to run from the Matlab command line for a standard mediation analysis.  

The toolbox also has special functions for Mediation Effect Parametric Mapping, the practice of running mediation on each voxel in a neuroimaging dataset (single or multilevel) and saving maps of mediation effects:
mediation_brain                 % for single-level mediation
mediation_brain_multilevel      % for multilevel mediation

The toolbox also has the ability to perform Multivariate mediation to identify a pattern across dense, high-dimensional mediators. See the bibliography below, and other tutorials in this series, for examples of voxel-wise mediation effect mapping and multivariate mediation.

**Installation and dependencies**

- Clone or download the [Mediation toolbox](https://github.com/canlab/MediationToolbox)
- Clone or download the [CANlab Core Tools toolbox](https://github.com/canlab/CanlabCore)
- Install SPM12 software (needed for neuroimaging mediation analyses only)
- Requires the Matlab statistics toolbox
- Add the toolboxes above "with subfolders" on your Matlab path
- See [canlab.github.io](https://canlab.github.io) for more detailed help installing CANlab Matlab toolboxes

**Tutorials**

Tutorials are included in the Mediation_walkthrough subfolder in the [Mediation Toolbox](https://github.com/canlab/MediationToolbox), and in the mediation analysis section of [https://canlab.github.io/walkthroughs](https://canlab.github.io/walkthroughs/)

Example datasets for mediation brain walkthroughs are included in the [Mediation Toolbox](https://github.com/canlab/MediationToolbox) and [CANlab Core toolbox](https://github.com/canlab/CanlabCore)

Tutorials include:
1. **mediation_1_basics.mlx** : A live script describing the fundamentals of mediation analysis and how to do single-level mediation analyses on simulated data
2. **mediation_example_script_1.m** : A script describing how to run single-level Mediation Effect Parametric Mapping with a sample fMRI dataset
3. **mediation_example_script_2.m** : Another version of above, with expanded descriptions 
4. **mediation_brain_multilevel_walkthrough1.mlx** : A live script demonstrating multilevel mediation with voxel-wise search using a sample fMRI dataset

The [walkthroughs](https://canlab.github.io/walkthroughs/) also include tutorials on CANlab object-oriented analysis and visualization, which can be used to visualize and further analyze mediation toolbox output.

For more tutorials on fMRI, see also:
[https://canlab.github.io/tutorials](https://canlab.github.io/tutorials/)


**Acknowledgements**

This toolbox was developed with the generous support of the U.S. National Science Foundation (NSF 0631637, Multilevel mediation techniques for fMRI) to Tor Wager and Martin Lindquist.  We are grateful to Prof. Niall Bolger and Prof. Michael Sobel for helpful discussions and input.		

**Key References**		

The key papers describing the toolbox are below. It would be helpful to cite these papers when using the M3 toolbox. The manuscripts and supplementary information contain fairly complete descriptions of the model and statistical procedures used.

    Shrout PE, Bolger N (2002) Mediation in experimental and nonexperimental studies: New procedures and recommendations. Psychol Methods 7:422–445.
 _   A classic reference for bootstrap-based inference in mediation._
    
    Kenny DA, Korchmaros JD, Bolger N (2003) Lower level mediation in multilevel models. Psychol Methods 8:115–128.
_    A seminal paper on multi-level mediation._
    
    Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., & Ochsner, K. N. (2008). Prefrontal-subcortical pathways mediating successful emotion regulation. Neuron, 59(6), 1037-1050. doi:10.1016/j.neuron.2008.09.006

_    This paper describes Mediation Effect Parametric Mapping with bootstrap-based inference and examples using multiple mediators. It describes and provides the first application of the canonical case where M is a set of brain images, and X and Y are given. It also describes suppressor effects and search when the initial variable (X) is the brain variable, and M and Y are given._
    
    Wager, T. D., Waugh, C. E., Lindquist, M., Noll, D. C., Fredrickson, B. L., & Taylor, S. F. (2009). Brain mediators of cardiovascular responses to social threat, Part I: Reciprocal dorsal and ventral sub-regions of the medial prefrontal cortex and heart-rate reactivity. NeuroImage, 47, 821-835.
 _    This paper describes multilevel mediation as applied to fMRI time series data, and provides the first application of multilevel mediation to fMRI data._
     
    Wager, T. D., van Ast, V. A., Hughes, B. L., Davidson, M. L., Lindquist, M. A., & Ochsner, K. N. (2009). Brain mediators of cardiovascular responses to social threat, Part II: Prefrontal-subcortical pathways and relationship with anxiety. NeuroImage, 47, 836-851.

_This paper describes multi-level mediation on fMRI time series, including the first application of moderated mediation and analyses of time-lagged mediators (mediators for which the time constants in relation to physiological outcomes differ across brain regions)._

    Atlas, L. Y., Bolger, N., Lindquist, M. A., & Wager, T. D. (2010). Brain mediators of predictive cue effects on perceived pain. J Neurosci, 30(39), 12964-12977. doi:10.1523/JNEUROSCI.0057-10.2010
   
_This paper applies multi-level mediation to single-trial brain images in fMRI, connecting an experimental design variable (X), brain mediators (M), and a behavioral outcome (Y), including second-level moderators. It provides the first application of multilevel mediation to single-trial fMRI data, and includes an extensive supplementary information document with more details on multilevel mediation._
   
    Woo, C. W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS biology, 13(1), e1002036. doi:10.1371/journal.pbio.1002036

_ This paper applies multi-level mediation to single-trial fMRI data, and provides the first application of multi-path (3 path) mediation, developed by Choong-Wan Woo._
        
    Oliver Y Chén, Ciprian Crainiceanu, Elizabeth L Ogburn, Brian S Caffo, Tor D Wager, Martin A Lindquist, High-dimensional multivariate mediation with application to neuroimaging data, Biostatistics, Volume 19, Issue 2, April 2018, Pages 121–136.

_    This paper describes the statistical foundations of multivariate mediation, and is the primary statistical reference for this technique._

    Stephan Geuter, Elizabeth A Reynolds Losin, Mathieu Roy, Lauren Y Atlas, Liane Schmidt, Anjali Krishnan, Leonie Koban, Tor D Wager, Martin A Lindquist, Multiple Brain Networks Mediating Stimulus–Pain Relationships in Humans, Cerebral Cortex, Volume 30, Issue 7, July 2020, Pages 4204–4219
    
_    This paper describes the application of multivariate mediation to a large single-trial fMRI dataset, and illustrates how mediation can be combined with other tools to interpret high-dimensional patterns of brain mediators._

: mediationanalysis.com, mediationanalysis.org

=======

**Additional References**

Here is a sample of additional papers using the M3 mediation toolbox. If you know of others we should include, please let us know!

Kober, H., Barrett, L. F., Joseph, J., Bliss-Moreau, E., Lindquist, K., & Wager, T. D. (2008). Functional grouping and cortical-subcortical interactions in emotion: A meta-analysis of neuroimaging studies. NeuroImage, 42, 998-1031.
 
Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., & Ochsner, K. N. (2008). Prefrontal-subcortical pathways mediating successful emotion regulation. Neuron, 59(6), 1037-1050. doi:10.1016/j.neuron.2008.09.006

Lim, S.-L., Padmala, S., & Pessoa, L. (2009). Segregating the significant from the mundane on a moment-to-moment basis via direct and indirect amygdala contributions. Proceedings of the National Academy of Sciences, 106(39), 16841-16846. 

Wager, T. D., Waugh, C. E., Lindquist, M., Noll, D. C., Fredrickson, B. L., & Taylor, S. F. (2009). Brain mediators of cardiovascular responses to social threat, Part I: Reciprocal dorsal and ventral sub-regions of the medial prefrontal cortex and heart-rate reactivity. NeuroImage, 47, 821-835. 

Wager, T. D., van Ast, V. A., Hughes, B. L., Davidson, M. L., Lindquist, M. A., & Ochsner, K. N. (2009). Brain mediators of cardiovascular responses to social threat, Part II: Prefrontal-subcortical pathways and relationship with anxiety. NeuroImage, 47, 836-851. 

Atlas, L. Y., Bolger, N., Lindquist, M. A., & Wager, T. D. (2010). Brain mediators of predictive cue effects on perceived pain. J Neurosci, 30(39), 12964-12977. doi:10.1523/JNEUROSCI.0057-10.2010

Buhle, J., & Wager, T. D. (2010). Performance-dependent inhibition of pain by an executive working memory task. Pain, 149(1), 19-26. doi:10.1016/j.pain.2009.10.027

Kober, H., Mende-Siedlecki, P., Kross, E. F., Weber, J., Mischel, W., Hart, C. L., & Ochsner, K. N. (2010). Prefrontal–striatal pathway underlies cognitive regulation of craving. Proceedings of the National Academy of Sciences, 107(33), 14811-14816. 

Krawitz, A., Fukunaga, R., & Brown, J. W. (2010). Anterior insula activity predicts the influence of positively framed messages on decision making. Cognitive, Affective, & Behavioral Neuroscience, 10(3), 392-405. 

Chen, C., Yang, C.-Y., & Cheng, Y. (2012). Sensorimotor resonance is an outcome but not a platform to anticipating harm to others. Social neuroscience, 7(6), 578-590.
 
Johnston, N. E., Atlas, L. Y., & Wager, T. D. (2012). Opposing effects of expectancy and somatic focus on pain. PLoS One, 7(6), e38854. doi:10.1371/journal.pone.0038854

Schneider, S., Peters, J., Bromberg, U., Brassen, S., Miedl, S. F., Banaschewski, T., . . . Garavan, H. (2012). Risk taking and the adolescent reward system: a potential common link to substance abuse. American Journal of Psychiatry. 

Atlas, L. Y., Lindquist, M. A., Bolger, N., & Wager, T. D. (2014). Brain mediators of the effects of noxious heat on pain. Pain, 155(8), 1632-1648. doi:10.1016/j.pain.2014.05.015

Denny, B. T., Ochsner, K. N., Weber, J., & Wager, T. D. (2014). Anticipatory brain activity predicts the success or failure of subsequent emotion regulation. Social cognitive and affective neuroscience, 9(4), 403-411. doi:10.1093/scan/nss148

Fan, Y.-T., Chen, C., Chen, S.-C., Decety, J., & Cheng, Y. (2014). Empathic arousal and social understanding in individuals with autism: evidence from fMRI and ERP measurements. Social cognitive and affective neuroscience, 9(8), 1203-1213. 

Cremers, H. R., Veer, I. M., Spinhoven, P., Rombouts, S. A., Yarkoni, T., Wager, T. D., & Roelofs, K. (2015). Altered cortical-amygdala coupling in social anxiety disorder during the anticipation of giving a public speech. Psychol Med, 45(7), 1521-1529. doi:10.1017/S0033291714002657

Woo, C. W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS biology, 13(1), e1002036. doi:10.1371/journal.pbio.1002036

Yamamoto, D. J., Woo, C. W., Wager, T. D., Regner, M. F., & Tanabe, J. (2015). Influence of dorsolateral prefrontal cortex and ventral striatum on risk avoidance in addiction: a mediation analysis. Drug Alcohol Depend, 149, 10-17. doi:10.1016/j.drugalcdep.2014.12.026

López-Solà, M., Woo, C.-W., Pujol, J., Deus, J., Harrison, B. J., Monfort, J., & Wager, T. D. (2016). Towards a neurophysiological signature for fibromyalgia. Pain. 

Vachon-Presseau, E., Roy, M., Woo, C.-W., Kunz, M., Martel, M.-O., Sullivan, M. J., . . . Rainville, P. (2016). Multiple faces of pain: Effects of chronic pain on the brain regulation of facial expression. Pain. 

van Ast, V. A., Spicer, J., Smith, E. E., Schmer-Galunder, S., Liberzon, I., Abelson, J. L., & Wager, T. D. (2016). Brain Mechanisms of Social Threat Effects on Working Memory. Cerebral cortex, 26(2), 544-556. doi:10.1093/cercor/bhu206

