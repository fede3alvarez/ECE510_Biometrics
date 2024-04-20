Overview:

	This porject uses Convulational Neural Networks in Matlab to perform Face Recognition.


Face Database:
	
	CalTech: https://data.caltech.edu/records/6rjah-hdv18
	Yale:	 https://www.kaggle.com/datasets/kpvisionlab/tufts-face-database

	During the development of this project, multiple Face Databases were visited and considered.
	In the end, the faces dataset used is a curated version of the Tufts and CalTech databases above.
	Images have been inspected to ensure that each subject has at least 10 images (Tufts) - through most of them are on the range of 20 images (CalTech).
	Any subject or images that did not meet this criteria was removed

	Additional, each images for a particular subject with separated in folders with arbitrary names to be used as labels / classes in the neural network training.
	

	Other Consideration

	https://github.com/microsoft/DigiFace1M
	Discarded due to the unsettling apperance of the subjects.
	Additionally, this dataset involved a lot more clean-up than the other datasets identified.
	
	https://www.kaggle.com/datasets/olgabelitskaya/yale-face-database
	Discarded due to the low number of images per subject.
	

Paper Considerations:
	
	Multiple papers where considered in the creation of this project.
	They can be found in zip files submitted with the deliverables.
	Multiple websites, and youtube channels were used when it comes to issues with code and syntax (those are not referenced here).
