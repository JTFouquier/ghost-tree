#!/usr/bin/env python
import click
from skbio.stats.distance import mantel
from skbio import TreeNode


@click.command()
@click.argument("tree_fp1")
@click.argument("tree_fp2")
@click.option("--method", default="pearson", type=click.Choice(["spearman", "pearson"]), help="correlation method to use in Mantel test")
def compare_trees(tree_fp1, tree_fp2, method):

    itstree = TreeNode.read(tree_fp1)
    hybridtree = TreeNode.read(tree_fp2)

    itstreedm = itstree.tip_tip_distances()
    hybridtreedm = hybridtree.tip_tip_distances()


    coeff, p_value, n = mantel(itstreedm,hybridtreedm,strict=False, method=method)

    click.echo("Correlation coefficient: %f" % coeff)
    click.echo("P-value: %f" % p_value)
    click.echo("Number of overlapping tips: %d" % n)


if __name__ == "__main__":
    compare_trees()
