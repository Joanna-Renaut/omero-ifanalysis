import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from itertools import combinations
from pathlib import Path
from ifanalysis.normalisation import normalise
from ifanalysis.sl_analysis import combplot_hist, intensity_plot
from ifanalysis._helper_functions import save_fig

STYLE = Path('~/matplotlib_style/Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


def custom_sort(val):
    order = {'NT': 0, 'DMSO': 1, 'iMAX only': 2, 'Pool': 3, 'RAD51': 4, 'BCL2': 5}
    return order.get(val, 6), val


def count_per_cond(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function to generate counts per condition and cell line data using groupby
    of the single cell dataframe from omero-screen. Data are grouped by
    cell line and condition. The function also normalises
    the data using normalise count function with the supplied ctr_cond as a reference.
    :param df: dataframe from omero-screen
    :param ctr_cond: control condition
    :return: dataframe with counts per condition and cell line
    """

    return (
        df.groupby(["experiment", "cell_line", "condition", "well", "plate_id", "genotype"])
        .size()
        .reset_index(name="abs_cell_count")
    )


if __name__ == '__main__':
    base_dir = Path.cwd()

    bio_repeat_folders = ['repeat-1', 'repeat-2', 'repeat-3']

    final_data_dir = base_dir / 'final-data'
    final_data_dir.mkdir(exist_ok=True)

    for repeat_folder in bio_repeat_folders:
        df = pd.DataFrame()
        repeat_path = base_dir / repeat_folder
        repeat_number = repeat_folder.split('-')[-1]

        EXPERIMENT_NAME = f'bio_repeat_0{repeat_number}'

        CONTROL_GENOTYPE = 'WT'
        SELECTED_CONDITIONS = []
        ALL_CONTROLS = ['DMSO']  # enter your controls here in the order you would like them in figures
        INTENSITY_PLOT_COLUMNS = []  # enter your column names for the different stains for an intensity plot
        RATIO_THRESHOLD = 0.7  # set the threshold for ratios shown on plot for relative cell counts

        for file in repeat_path.glob('*.csv'):
            data = pd.read_csv(file, index_col=0)
            data['experiment'] = EXPERIMENT_NAME
            df = pd.concat([df, data])
        df.reset_index(drop=True, inplace=True)

        df_counts = count_per_cond(df)

        # Process each experiment only once
        unique_experiments = df_counts['experiment'].unique()

        processed_conditions = set()

        # Filter out control and experimental groups
        control_df = df_counts[df_counts['genotype'] == CONTROL_GENOTYPE]

        grouped = df_counts.groupby(['experiment', 'condition', 'genotype'])['abs_cell_count']
        df_counts['abs_cell_count_sem'] = grouped.transform(lambda x: x.sem())

        # Merge experimental data with control data, including WT with itself
        merged_df = pd.merge(df_counts, control_df[['experiment', 'condition', 'abs_cell_count']],
                             on=['experiment', 'condition'], suffixes=('', '_ctrl'))

        # Calculate ratios, setting WT against itself to have a meaningful value (e.g., 1)
        merged_df['ratio'] = np.where(merged_df['genotype'] == CONTROL_GENOTYPE, 1,
                                      merged_df['abs_cell_count'] / merged_df['abs_cell_count_ctrl'])

        # Prepare the final DataFrame for plotting
        final_plot_df = merged_df[['experiment', 'condition', 'genotype', 'ratio', 'abs_cell_count']].copy()
        final_plot_df.rename(columns={'genotype_exp': 'genotype'}, inplace=True)

        # Sorting and selection based on conditions
        sorted_values = sorted(merged_df['condition'].unique(), key=custom_sort)
        SELECTED_CONDITIONS = sorted_values

        # Sort by condition using custom_sort for plotting
        final_plot_df['condition'] = pd.Categorical(final_plot_df['condition'], categories=sorted_values, ordered=True)
        final_plot_df.sort_values(by=['condition', 'genotype'], inplace=True)

        # Print out the first two colors for verification
        print("Color at index 0:", colors[0])
        print("Color at index 1:", colors[1])
        unique_genotypes = final_plot_df['genotype'].unique()

        # Check order of categories in your data
        print("Unique genotypes in order:", unique_genotypes)

        selected_conds = SELECTED_CONDITIONS

        if INTENSITY_PLOT_COLUMNS:
            columns = INTENSITY_PLOT_COLUMNS

            intensity_plot(df, selected_conds, columns, title='Intensity Plot, WT', hue='genotype')

        # Get unique genotypes
        unique_genotypes = list(final_plot_df['genotype'].unique())

        if CONTROL_GENOTYPE in unique_genotypes:
            unique_genotypes.remove(CONTROL_GENOTYPE)
        hue_order = [CONTROL_GENOTYPE] + unique_genotypes

        # create figures for each genotype (cell line) in the experiment
        for genotype in unique_genotypes:

            df1 = df[df.genotype == 'WT']
            df2 = df[df.genotype == genotype]

            df1_cc = normalise(df1, ['integrated_int_DAPI'])
            df1_cc['integrated_int_DAPI_norm'] = df1_cc['integrated_int_DAPI_norm'] * 2
            df2_cc = normalise(df2, ['integrated_int_DAPI'])
            df2_cc['integrated_int_DAPI_norm'] = df2_cc['integrated_int_DAPI_norm'] * 2

            combplot_hist(df1_cc, selected_conds, 'CellCycle Histogram, WT')
            combplot_hist(df2_cc, selected_conds, f'CellCycle Histogram, {genotype}')

            sorted_values = sorted(final_plot_df.condition.unique(), key=custom_sort)

            fig, ax = plt.subplots(figsize=(6, 6))
            sns.barplot(x='condition', y='abs_cell_count', hue='genotype', data=final_plot_df, ax=ax,
                        hue_order=hue_order, order=sorted_values)
            plt.xticks(rotation=90)
            ax.legend(loc='upper right')

            bar_width = 0.8

            correction_factor = bar_width / 2.0 * (len(hue_order) - 1)

            for i, geno in enumerate(hue_order):
                subset = final_plot_df[final_plot_df['genotype'] == geno]
                positions = np.arange(len(sorted_values)) + (i * bar_width) - correction_factor

                # Check if positions and subset lengths match before plotting
                if not subset.empty and len(positions) == len(subset):
                    # Correctly calculate positions for this genotype's bars
                    positions = np.arange(len(sorted_values)) + (
                        i - len(hue_order) / 2.0) * 0.2  # Adjust the 0.2 if necessary to align with your bar width
                    # Plot error bars
                    plt.errorbar(positions, subset['abs_cell_count'], yerr=subset['abs_cell_count_sem'], fmt='none',
                                 color='black', capsize=3, elinewidth=1)

            save_fig(fig, repeat_path, f'screen{repeat_number}_hits_abscounts')
            plt.savefig(f'{repeat_path}/screen{repeat_number}_hits_abscounts.png', dpi=600)

            # Sort by condition using custom_sort for plotting
            sorted_conditions = sorted(final_plot_df['condition'].unique(), key=custom_sort)
            final_plot_df['condition'] = pd.Categorical(final_plot_df['condition'], categories=sorted_conditions,
                                                        ordered=True)
            final_plot_df.sort_values(by=['condition', 'genotype'], inplace=True)

            final_genotype = final_plot_df[final_plot_df['genotype'] == genotype].copy()

            bar_colors = [colors[1] if ratio < RATIO_THRESHOLD else colors[0] for ratio in
                          final_genotype.groupby('condition', observed=True)['ratio'].mean()]

            # Create bar plot for ratios
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.barplot(x='condition', y='ratio', data=final_genotype, ax=ax, palette=bar_colors,
                        order=sorted_conditions)

            plt.axhline(y=RATIO_THRESHOLD, color='gray', linestyle='--', linewidth=2)
            plt.xticks(rotation=90)
            plt.tight_layout()

            # Save the figure
            fig_filename = f'screen{repeat_number}_ratios_{genotype}'
            save_fig(fig, repeat_path, fig_filename)
            plt.savefig(f'{repeat_path}/{fig_filename}.png', dpi=600)
            plt.close(fig)

        # Aggregate abs_cell_count by experiment, condition, and genotype to get mean values
        aggregated_df = df_counts.groupby(['experiment', 'condition', 'genotype'], as_index=False).agg({
            'abs_cell_count': 'mean'
        })

        # Merge aggregated data with control data to calculate ratios
        control_counts = aggregated_df[aggregated_df['genotype'] == CONTROL_GENOTYPE][
            ['experiment', 'condition', 'abs_cell_count']].rename(columns={'abs_cell_count': 'control_count'})
        merged_aggregated_df = pd.merge(aggregated_df, control_counts, on=['experiment', 'condition'])

        # Calculate ratio for non-control genotypes, set ratio as 1 for control genotype
        merged_aggregated_df['ratio'] = np.where(merged_aggregated_df['genotype'] == CONTROL_GENOTYPE, 1,
                                                 merged_aggregated_df['abs_cell_count'] / merged_aggregated_df[
                                                     'control_count'])

        # Now, merged_aggregated_df contains aggregated abs_cell_count and calculated ratios
        final_data = merged_aggregated_df[['experiment', 'condition', 'genotype', 'abs_cell_count', 'ratio']].copy()

        # Sorting by condition as per custom sort for plotting and analysis purposes
        final_data['condition'] = pd.Categorical(final_data['condition'], categories=sorted_values, ordered=True)
        final_data.sort_values(by=['condition', 'genotype'], inplace=True)

        final_data.to_csv(final_data_dir / f'bio_repeat_0{repeat_number}_analysis_results.csv', index=False)
        print(f"Results saved to 'bio_repeat_0{repeat_number}_analysis_results.csv'.")

    stats_dir = base_dir / 'stats'
    stats_dir.mkdir(exist_ok=True)
    final_files = final_data_dir.glob('*.csv')
    analyses_data = pd.concat([pd.read_csv(f) for f in final_files])
    analyses_data['abs_cell_count'] = pd.to_numeric(analyses_data['abs_cell_count'], errors='coerce')
    analyses_data['ratio'] = pd.to_numeric(analyses_data['ratio'], errors='coerce')

    # Aggregate data to calculate means
    cell_count_means = analyses_data.groupby(['experiment', 'condition'])['abs_cell_count'].mean().reset_index()
    cell_count_means.to_csv(stats_dir / 'means_of_abs_cell_count.csv', index=False)
    ratio_means = analyses_data.groupby(['experiment', 'condition'])['ratio'].mean().reset_index()
    ratio_means.to_csv(stats_dir / 'means_of_relative_cell_count.csv', index=False)

    # Calculate standard deviations
    cell_count_stds = analyses_data.groupby(['experiment', 'condition'])['abs_cell_count'].std().reset_index()
    cell_count_stds.to_csv(stats_dir / 'std_of_abs_cell_count.csv', index=False)
    ratio_stds = analyses_data.groupby(['experiment', 'condition'])['ratio'].std().reset_index()
    ratio_stds.to_csv(stats_dir / 'std_of_relative_cell_count.csv', index=False)
    print('Means and std saved to csv')

    # Initialize a list to store t-test results
    t_test_results = []
    # Iterate over each experiment
    for experiment in analyses_data['experiment'].unique():
        experiment_data = analyses_data[analyses_data['experiment'] == experiment]
        # Identify all unique conditions within this experiment
        conditions = experiment_data['condition'].unique()

        # Generate all possible pairs of conditions
        for condition1, condition2 in combinations(conditions, 2):
            # Extract the abs_cell_count for each condition
            data1 = experiment_data[experiment_data['condition'] == condition1]['abs_cell_count']
            data2 = experiment_data[experiment_data['condition'] == condition2]['abs_cell_count']

            # Perform t-test between the two conditions
            t_stat, p_val = ttest_ind(data1, data2, equal_var=False)  # assuming unequal variance

            # Store the result
            t_test_results.append({
                'experiment': experiment,
                'condition_pair': f"{condition1} vs {condition2}",
                'p_value': p_val
            })

    # Convert the list of results into a DataFrame
    p_values_df = pd.DataFrame(t_test_results)
    print(p_values_df)
    p_values_df.to_csv(stats_dir / 'p-values.csv', index=False)
    print('p-values saved to csv. Analysis complete.')
