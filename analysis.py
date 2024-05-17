import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from pathlib import Path
STYLE = Path('Style_01.mplstyle')
plt.style.use(STYLE)
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


def calculate_sem(group):
    return np.std(group) / np.sqrt(len(group))


def add_significance_markers(ax, p_values, x_positions, data, offset=0.05):
    for x, p in zip(x_positions, p_values):
        if p < 0.001:
            ax.text(x, data.iloc[x]['average_ratio'] + data.iloc[x]['sem'] + offset, '***', ha='center', va='bottom',
                    color='black', fontsize=15)
        elif p < 0.01:
            ax.text(x, data.iloc[x]['average_ratio'] + data.iloc[x]['sem'] + offset, '**', ha='center', va='bottom',
                    color='black', fontsize=15)
        elif p < 0.05:
            ax.text(x, data.iloc[x]['average_ratio'] + data.iloc[x]['sem'] + offset, '*', ha='center', va='bottom',
                    color='black', fontsize=15)


def plot_ratios_with_error_bars(data, ratio_threshold, control_genotype, ko_genotype, all_controls, combined_data):
    # Perform Welch's t-tests and get p-values
    p_values = []
    for condition in data['condition'].unique():
        ko_values = combined_data[(combined_data['condition'] == condition) &
                                  (combined_data['genotype'] == ko_genotype)]['normalized_count']
        wt_values = combined_data[(combined_data['condition'] == condition) &
                                  (combined_data['genotype'] == control_genotype)]['normalized_count']

        # # Add a small amount of noise to avoid precision loss (not recommended but can bypass the warning)
        # ko_values += np.random.normal(0, 1e-8, ko_values.shape)
        # wt_values += np.random.normal(0, 1e-8, wt_values.shape)
        print(f"Condition: {condition}")  # Debugging: Print condition
        print(f"KO values: {ko_values.tolist()}")  # Debugging: Print KO values
        print(f"WT values: {wt_values.tolist()}")  # Debugging: Print WT values

        if len(ko_values) > 1 and len(wt_values) > 1:  # Ensure there are enough data points
            _, p_value = ttest_ind(ko_values, wt_values, equal_var=False)
            if np.isnan(p_value):  # Check if precision loss occurred
                ko_values += np.random.normal(0, 1e-8, ko_values.shape)
                wt_values += np.random.normal(0, 1e-8, wt_values.shape)
                _, p_value = ttest_ind(ko_values, wt_values, equal_var=False)
            p_values.append(p_value)
        else:
            p_values.append(np.nan)  # Append NaN if there's not enough data

    # test p-values for code testing
    # p_values = [0.0001, 0.005, 0.03, 0.2, 0.04, 0.02, 0.0005, 0.07, 0.15, 0.001]

    # Sort conditions: controls first, then others alphabetically
    sorted_conditions = all_controls + sorted([cond for cond in data['condition'].unique() if cond not in all_controls])

    # Reorder data based on sorted conditions
    data = data.set_index('condition').loc[sorted_conditions].reset_index()

    # Print data and p-values for debugging
    # print(data.head(70))
    print(p_values)

    # Plot
    plt.figure(figsize=(12, 8))
    # colors = ['red' if ratio < ratio_threshold else 'blue' for ratio in data['average_ratio']]
    bar_color = [colors[1] if ratio < RATIO_THRESHOLD else colors[0] for ratio in data['average_ratio']]
    bars = plt.bar(data['condition'], data['average_ratio'], yerr=data['sem'], color=bar_color, capsize=5,
                   edgecolor='black')

    # Set the order of conditions and rotate labels
    plt.xticks(ticks=np.arange(len(sorted_conditions)), labels=sorted_conditions, rotation=90, ha='center')

    # Add dashed line for ratio threshold
    plt.axhline(y=ratio_threshold, color='gray', linestyle='--')

    # Add significance markers
    add_significance_markers(plt.gca(), p_values, range(len(data)), data, offset=0.03)

    plt.xlabel('Condition')
    plt.ylabel('Averaged KO/WT Ratio')
    plt.title(f'Normalized KO/WT Ratios for {ko_genotype}')

    handles = [
        plt.Line2D([0], [0], color='none', marker='o', linestyle='none', markersize=10, label='* p < 0.05'),
        plt.Line2D([0], [0], color='none', marker='o', linestyle='none', markersize=10, label='** p < 0.01'),
        plt.Line2D([0], [0], color='none', marker='o', linestyle='none', markersize=10, label='*** p < 0.001')
    ]
    plt.legend(handles=handles, loc='upper right', title='Significance', handletextpad=-1, labelspacing=1.2)

    plt.tight_layout()

    figure_path = base_dir / 'figures'
    figure_path.mkdir(exist_ok=True)
    plt.savefig(f'{figure_path}/ratio_plot_for_{ko_genotype}_all_repeats.png', dpi=600, format='png')
    plt.savefig(f'{figure_path}/ratio_plot_for_{ko_genotype}_all_repeats.pdf', format='pdf')
    plt.show()

    # Save p-values to CSV
    p_values_df = pd.DataFrame({'condition': sorted_conditions, 'p-value': p_values})
    p_values_df.to_csv(final_data_dir / f'p_values_{ko_genotype}.csv', index=False)


if __name__ == '__main__':
    base_dir = Path.cwd()

    bio_repeat_folders = ['repeat-1', 'repeat-2', 'repeat-3']

    final_data_dir = base_dir / 'final-data'
    final_data_dir.mkdir(exist_ok=True)

    all_data = []
    CONTROL_GENOTYPE = 'WT'
    KO_GENOTYPES = ['GENE1', 'GENE2']  # List of KO genotypes
    MAIN_CONTROL = 'NT'
    ALL_CONTROLS = ['NT', 'DMSO']  # enter your controls here in the order you would like them in figures
    INTENSITY_PLOT_COLUMNS = []  # enter your column names for the different stains for an intensity plot
    RATIO_THRESHOLD = 0.7  # set the threshold for ratios shown on plot for relative cell counts

    for repeat_folder in bio_repeat_folders:
        df = pd.DataFrame()
        repeat_path = base_dir / repeat_folder
        repeat_number = repeat_folder.split('-')[-1]

        EXPERIMENT_NAME = f'bio_repeat_0{repeat_number}'

        for file in repeat_path.glob('*.csv'):
            data = pd.read_csv(file, index_col=0, low_memory=False)
            data['experiment'] = EXPERIMENT_NAME
            df = pd.concat([df, data])
        df.reset_index(drop=True, inplace=True)
        df['condition'] = df['condition'].astype(str)

        # Calculate the cell counts per well
        cell_counts = df.groupby(
            ['experiment', 'plate_id', 'well', 'cell_line', 'condition', 'genotype']).size().reset_index(
            name='cell_count')

        # Average technical replicates
        avg_cell_counts = cell_counts.groupby(
            ['experiment', 'cell_line', 'condition', 'genotype']).cell_count.mean().reset_index(name='tech_repeat_mean')

        # Calculate SEM for error bars
        avg_cell_counts['sem'] = cell_counts.groupby(['experiment', 'cell_line', 'condition', 'genotype'])[
            'cell_count'].transform(lambda x: x.sem())

        # Calculate the mean values for NT condition
        nt_mean_ko = avg_cell_counts[(avg_cell_counts['genotype'] == KO_GENOTYPES[0]) &
                                     (avg_cell_counts['condition'] == MAIN_CONTROL)]['tech_repeat_mean'].mean()
        nt_mean_wt = avg_cell_counts[(avg_cell_counts['genotype'] == CONTROL_GENOTYPE) &
                                     (avg_cell_counts['condition'] == MAIN_CONTROL)]['tech_repeat_mean'].mean()

        # Normalize the treated counts by the NT control counts
        avg_cell_counts['normalized_count'] = avg_cell_counts.apply(
            lambda row: row['tech_repeat_mean'] / nt_mean_ko if row['genotype'] in KO_GENOTYPES else
            row['tech_repeat_mean'] / nt_mean_wt, axis=1)

        # Normalize the SEM values
        avg_cell_counts['normalized_sem'] = avg_cell_counts.apply(
            lambda row: row['sem'] / nt_mean_ko if row['genotype'] in KO_GENOTYPES else
            row['sem'] / nt_mean_wt, axis=1)

        all_data.append(avg_cell_counts)

    # Combine normalized data from all repeats
    combined_data = pd.concat(all_data)

    for ko_genotype in KO_GENOTYPES:
        # Group by condition and genotype to calculate the average ratio and SEM across repeats
        final_data = combined_data[combined_data['genotype'].isin([
            CONTROL_GENOTYPE, ko_genotype])].groupby(['condition', 'genotype']).agg(
            {'normalized_count': 'mean', 'normalized_sem': 'mean'}).reset_index()

        # Pivot to get KO and WT in separate columns for easier ratio calculation
        final_data_pivot = final_data.pivot(index='condition', columns='genotype',
                                            values=['normalized_count', 'normalized_sem']).reset_index()
        final_data_pivot.columns = ['condition', f'normalized_count_{ko_genotype}', 'normalized_count_WT',
                                    f'sem_{ko_genotype}', 'sem_WT']

        # Calculate average ratio and SEM
        final_data_pivot['average_ratio'] = (final_data_pivot[f'normalized_count_{ko_genotype}'] /
                                             final_data_pivot['normalized_count_WT'])
        final_data_pivot['sem'] = final_data_pivot[[f'sem_{ko_genotype}', 'sem_WT']].max(axis=1)  # Combine SEMs

        # Reorder data based on sorted conditions
        sorted_conditions = ALL_CONTROLS + sorted([cond for cond in final_data_pivot['condition'].unique() if
                                                   cond not in ALL_CONTROLS])
        final_data_pivot = final_data_pivot.set_index('condition').loc[sorted_conditions].reset_index()

        # Save final_data to CSV
        final_data_pivot.to_csv(final_data_dir / f'final_data_{ko_genotype}.csv', index=False)

        # Plot ratios with error bars and save p-values
        plot_ratios_with_error_bars(final_data_pivot, RATIO_THRESHOLD, CONTROL_GENOTYPE, ko_genotype,
                                    ALL_CONTROLS, combined_data)
