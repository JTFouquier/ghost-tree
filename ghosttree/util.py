from skbio import TreeNode
from skbio.stats.distance import mantel


def compare_tip_to_tip_distances(tree_fh1, tree_fh2, method="pearson"):
    tree1 = TreeNode.read(tree_fh1)
    tree2 = TreeNode.read(tree_fh2)

    dm1 = tree1.tip_tip_distances()
    dm2 = tree2.tip_tip_distances()

    return mantel(dm1, dm2, strict=False, method=method)
