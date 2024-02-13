import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from itertools import combinations
from pathlib import Path
from ifanalysis.counts import count_cells
from ifanalysis.normalisation import normalise
from ifanalysis.sl_analysis import combplot_hist, intensity_plot
from ifanalysis._helper_functions import save_fig
STYLE = Path('~/matplotlib_style/Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


def custom_sort(val):
    if val == 'NT':
        return (0, val)
    elif val == 'iMAX only':
        return (1, val)
    elif val == 'Pool':
        return (2, val)
    else:
        return (3, val)


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
        .reset_index(name="abs cell count")
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
        for file in repeat_path.glob('*.csv'):
            data = pd.read_csv(file, index_col=0)
            data['experiment'] = EXPERIMENT_NAME
            df = pd.concat([df, data])
        df.reset_index(drop=True, inplace=True)

        df_counts = count_per_cond(df)
        results = []

        # Process each experiment only once
        unique_experiments = df_counts['experiment'].unique()

        processed_conditions = set()

        for experiment in unique_experiments:

            experiment_data = df_counts[df_counts['experiment'] == experiment]
            unique_conditions = experiment_data['condition'].unique()

            for condition in unique_conditions:
                if (experiment, condition) in processed_conditions:
                    continue
                condition_data = experiment_data[experiment_data['condition'] == condition]
                unique_genotypes = condition_data['genotype'].unique()

                # Control count for the condition within the experiment
                control_count = condition_data[condition_data['genotype'] == CONTROL_GENOTYPE]['abs cell count'].sum()
                processed_conditions.add((experiment, condition))

                # Now, loop through each genotype for the specific condition
                for genotype in unique_genotypes:
                    genotype_data = condition_data[condition_data['genotype'] == genotype]
                    genotype_count = genotype_data['abs cell count'].sum()
                    ratio = genotype_count / control_count if control_count > 0 else np.nan

                    # Append results for final DataFrame
                    results.append({
                        'experiment': experiment,
                        'condition': condition,
                        'genotype': genotype,
                        'abs_cell_count': genotype_count,
                        'ratio': ratio if genotype != CONTROL_GENOTYPE else 1
                    })
                    conditions_df = pd.DataFrame(results)
                    conditions_df = conditions_df[conditions_df['condition'].isin(['DMSO', 'inhibitor'])]
                    # print(conditions_df.condition.unique())
                    sorted_values = sorted(conditions_df.condition.unique(), key=custom_sort)
                    # print(sorted_values)

                    # Sort conditions_df based on sorted_values
                    conditions_df['condition'] = conditions_df['condition'].astype("category")
                    conditions_df['condition'] = conditions_df['condition'].cat.set_categories(sorted_values,
                                                                                               ordered=True)
                    conditions_df.sort_values(by='condition', inplace=True)

                    # Define custom colors for bars based on sorted conditions_df
                    bar_colors = [colors[1] if ratio < 0.7 else colors[0] for ratio in
                                  conditions_df.groupby('condition', observed=True)['ratio'].mean()]

                    # create bar plot
                    fig, ax = plt.subplots(figsize=(10, 6))
                    sns.barplot(x='condition', y='ratio', data=conditions_df, ax=ax, palette=bar_colors)
                    # ax.legend_.remove()
                    plt.xticks(rotation=90)
                    plt.title('Ratios by Condition')
                    plt.tight_layout()
                    save_fig(fig, repeat_path, f'{genotype}_screen{repeat_number}_ratios')

                    mean_ratios = conditions_df.groupby('condition', observed=True)['ratio'].mean()

                    conditions = [condition for condition in mean_ratios.index if mean_ratios[condition] < 0.7]
                    items_to_add = ['DMSO', 'inhibitor']
                    conditions = items_to_add + conditions
                    # print(conditions)

                    df_hits_counts = df_counts[df_counts['condition'].isin(conditions)]
                    # Print out the first two colors for verification
                    print("Color at index 0:", colors[0])
                    print("Color at index 1:", colors[1])

                    # Check order of categories in your data
                    print("Unique genotypes in order:", df_hits_counts['genotype'].unique())

                    fig, ax = plt.subplots(figsize=(6, 6))
                    color_mapping = {'WT': colors[0], f'{genotype}': colors[1]}
                    hue_order = ['WT', f'{genotype}']
                    sns.barplot(x='condition', y='abs cell count', data=df_hits_counts, order=conditions,
                                palette=color_mapping,
                                hue='genotype', ax=ax, hue_order=hue_order)
                    plt.xticks(rotation=90)
                    save_fig(fig, repeat_path, f'{genotype}_screen{repeat_number}_hits_abscounts')

                    selected_conds = ['DMSO', 'inhibitor']

                    columns = ['intensity_mean_Casp3/7_nucleus', 'intensity_mean_p21_nucleus',
                               'intensity_mean_yH2A.X_nucleus']
                    intensity_plot(df, selected_conds, columns, title='Intensity Plot, WT', hue='genotype')

                    df1 = df[df.genotype == 'WT']
                    df2 = df[df.genotype == genotype]

                    df1_cc = normalise(df1, ['integrated_int_DAPI'])
                    df1_cc['integrated_int_DAPI_norm'] = df1_cc['integrated_int_DAPI_norm'] * 2
                    df2_cc = normalise(df2, ['integrated_int_DAPI'])
                    df2_cc['integrated_int_DAPI_norm'] = df2_cc['integrated_int_DAPI_norm'] * 2

                    combplot_hist(df1_cc, selected_conds, 'CellCycle Histogram, WT')
                    combplot_hist(df2_cc, selected_conds, f'CellCycle Histogram, {genotype}')

        final_df = pd.DataFrame(results)
        final_df.to_csv(final_data_dir / f'bio_repeat_0{repeat_number}_analysis_results.csv', index=False)
        print(f"Results saved to 'bio_repeat_0{repeat_number}_analysis_results.csv'.")
        # print(final_df)

    stats_dir = base_dir / 'stats'
    stats_dir.mkdir(exist_ok=True)
    final_files = final_data_dir.glob('*.csv')
    analyses_data = pd.concat([pd.read_csv(f) for f in final_files])
    # print(analyses_data.head())
    analyses_data['abs_cell_count'] = pd.to_numeric(analyses_data['abs_cell_count'], errors='coerce')
    analyses_data['ratio'] = pd.to_numeric(analyses_data['ratio'], errors='coerce')
    # Step 2: Aggregate data to calculate means
    cell_count_means = analyses_data.groupby(['experiment', 'condition'])['abs_cell_count'].mean().reset_index()
    # print(cell_count_means)
    cell_count_means.to_csv(stats_dir / 'means_of_abs_cell_count.csv', index=False)
    ratio_means = analyses_data.groupby(['experiment', 'condition'])['ratio'].mean().reset_index()
    # print(ratio_means)
    ratio_means.to_csv(stats_dir / 'means_of_relative_cell_count.csv', index=False)
    # Step 3: Calculate standard deviations
    cell_count_stds = analyses_data.groupby(['experiment', 'condition'])['abs_cell_count'].std().reset_index()
    # print(cell_count_stds)
    cell_count_stds.to_csv(stats_dir / 'std_of_abs_cell_count.csv', index=False)
    ratio_stds = analyses_data.groupby(['experiment', 'condition'])['ratio'].std().reset_index()
    # print(ratio_stds)
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
