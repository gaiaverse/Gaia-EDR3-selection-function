.. external-libs

#####################
External Libraries
#####################


Our code ships with static copies of two external libraries:

* `Eigen <https://eigen.tuxfamily.org>`_, distributed under the appropriate open-source license
* `JSL <https://github.com/JTFraser/JSL>`_, a custom library also written by me, but kept elsewhere for reusability

The static nature of the shipped libraries, and the fact that the accompanying makefile points directly to them should make the code safe against whatever conflicting versions of these libraries you may have (though I don't fancy guaranteeing that). 
