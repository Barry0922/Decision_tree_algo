# Decision_tree_algo
Build decision tree algo by using c code.

Normally, people use sklearn module to use decision tree algo.
Here, by using c code to reproduce decision tree algo.

--------------------------------------------------------------------------------------------------------------

By using greedy search to find best features to split training data into sub-data.
Best features means it has highest information gain. In other words, the sub-data has low entropy (high purity).
Using tree structutre to record all relevant infaomtyion and using this tree as rules to identitfy data.

--------------------------------------------------------------------------------------------------------------

The iris samples dataset is from scikit-learn. The most famous ML samples dataset around the word.

iris03.csv is the iris data which has been splitted by using "train_test_split " from "sklearn.model_selection". (random state =1, test_size = 0.3)
iris04.csv is the iris data which has been splitted by using "train_test_split " from "sklearn.model_selection". (random state =1, test_size = 0.4)
iris03 and iris04 are the visualization of tree from sklearn.

After running this program, we can use above iris pdf file to check whether our program is right or not.
