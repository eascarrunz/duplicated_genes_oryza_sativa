from matplotlib import pyplot as plt
import pandas as pd

if __name__ == "__main__":
    fig, axes = plt.subplots(3, 1, sharey=True, sharex=True)

    for s in range(3):
        df_hs = pd.read_csv(f"output/hstags{s}.tsv", sep='\t')
        df_ls = pd.read_csv(f"output/lstags{s}.tsv", sep='\t')
        max_tag_size = max(df_hs["size"].max(), df_ls["size"].max())
        axes[s].set_title(f"Spacer tolerance = {s}", y=0.7)
        axes[s].hist(x=df_ls["size"], label="Low Stringency", alpha=0.5, bins=range(2, max_tag_size), log=True)
        axes[s].hist(x=df_hs["size"], label="High Stringency", alpha=0.5, bins=range(2, max_tag_size), log=True)
        axes[s].set_ylabel("Log frequency")
        
    axes[0].legend()
    axes[2].set_xlabel("Number of genes in TAG")
    fig.tight_layout()

    fig.savefig("output/tag_histograms.png", dpi=200)
