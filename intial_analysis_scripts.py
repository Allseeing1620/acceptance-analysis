
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from aa_helpers import create_plot_with_background
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

def concat_csvs_with_unique_events(files):
    """Load and concatenate CSV files with globally unique event IDs"""
    dfs = []
    offset = 0

    for file in files:
        df = pd.read_csv(file)
        df['event'] = df['event'] + offset
        offset = df['event'].max() + 1      # Set offset for next file
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


files = sorted(glob.glob('data\\k_lambda_5x41_5000evt_*.mcpart_lambda.csv.zip'))

if len(files) == 0:
    print("Check the path and file name pattern.")
else:

    combined_df = concat_csvs_with_unique_events(files)
    combined_df.to_feather('combined_all_lambda_data.feather')

df = pd.read_feather("combined_all_lambda_data_18x275.feather")

df_not_dec = df[(df['lam_is_first'] == 1) & (df['prot_id'].isna() & df['neut_id'].isna())]

df_pram = df[df['lam_is_first'] == 1]
df_not_pram = df[df['lam_is_first'] == 0]

# proton + pi-

p_pi_minus_decays = df_pram[
   df_pram['prot_id'].notna() & 
   df_pram['pimin_id'].notna() &
   df_pram['neut_id'].isna() &    
   df_pram['pizero_id'].isna()    
].copy()\

# neutron + pi0 -> gamma gamma 
n_pi_zero_decays = df_pram[
    df_pram['prot_id'].isna() & 
    df_pram['pimin_id'].isna() &
    df_pram['neut_id'].notna() &    
    df_pram['pizero_id'].notna() &
    df_pram['gamtwo_id'].notna() &
    df_pram['gamone_id'].notna()
].copy()



def calculate_decayed():
    
    total = len(df)

    proton_only = df['prot_id'].notna() & df['neut_id'].isna()
    neutron_only = df['neut_id'].notna() & df['prot_id'].isna()
    not_decayed = df['prot_id'].isna() & df['neut_id'].isna()
    
    p_proton = (proton_only.sum() / total) * 100
    p_neutron = (neutron_only.sum() / total) * 100  
    p_not_decayed = (not_decayed.sum() / total) * 100

    print(f"\n decayed proton: {p_proton:.1f}%")
    print(f"\ndecayed neutron: {p_neutron:.1f}%")
    print(f"\n not decayed: {p_not_decayed:.1f}%")


def calculate_percentage(df_all, condition_mask):
    total = len(df_all)
    percentage = (condition_mask.sum() / total) * 100
    return percentage

def analyze_primary_vs_secondary_lambdas():
    print(f"Secondary lambdas:{calculate_percentage(df, df['lam_is_first'] == 0):.1f}%")
    

def plot_point(x_axis, y_axis,xlabel = "z [mm]" ,ylabel = "y [mm]"):

    fig, ax = create_plot_with_background()
    # Optional: overlay the reference points to verify alignment
    ax.plot(x_axis, y_axis, marker="o", linestyle="none", alpha=0.3)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)

    plt.tight_layout()
    plt.show()


def plot_primary_lambda_decay_z_distribution():
    """1D гистограмма Z-координаты распада первичных лямбд"""
    plt.figure(figsize=(12, 6))
    plt.hist(df_pram['lam_epz'], bins=80, alpha=0.7, edgecolor='black')
    plt.xlabel('Z decay vertex coordinate of Λ⁰ (mm)')
    plt.ylabel('Number of events')
    plt.title('Distribution of Z decay vertex coordinate for primary Λ⁰')
    plt.grid(True, alpha=0.3)
    
    # Limit X-axis to
    plt.xlim(0, 45000)
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(5000))
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1000))
    
    plt.grid(True, which='major', alpha=0.3)
    plt.grid(True, which='minor', alpha=0.1)
    
    plt.tight_layout()
    plt.show()


def plt_hist2d(x_axis, y_axis,bins=50, xlabel = "z [mm]" ,ylabel = "y [mm]",title = "Title is missing"):

    fig, ax = create_plot_with_background()

    # Сохраняем оригинальные пределы осей
    original_xlim = ax.get_xlim()
    original_ylim = ax.get_ylim()

    # Строим 2D гистограмму для первичных лямбд
    h = ax.hist2d(x_axis,y_axis, bins=bins,
                   cmap='viridis', norm=LogNorm())

    # Восстанавливаем пределы осей
    ax.set_xlim(original_xlim)
    ax.set_ylim(original_ylim)

    # Добавляем colorbar
    fig.colorbar(h[3], ax=ax, label='Counts', shrink=0.3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)

    ax.set_title(title)

    plt.tight_layout()
    plt.show()


# 2Д гистограмма на фоне детектора, где распадаются лямбды (for proton+pi-)

# проблемы с масштабом 


# взять только 20 евентов (столько, сколько хорошо видно разницу) 
# построить так: красная паочка как летят лямбды (точка старт - ф ниш), синяя палочка - как летят протоны, зеленая - пионы
def decay_trajectories():
    sample_df = p_pi_minus_decays.iloc[50:200].head(3)
    
    # Координаты для лямбд (красные)
    lam_start_z = sample_df['lam_vz'].values  # начало - рождение
    lam_start_x = sample_df['lam_vx'].values
    lam_end_z = sample_df['lam_epz'].values   # конец - распад  
    lam_end_x = sample_df['lam_epx'].values
    
    # Координаты для протонов (синие) - летят от точки распада лямбды
    prot_start_z = sample_df['lam_epz'].values  # начало - точка распада лямбды
    prot_start_x = sample_df['lam_epx'].values  
    prot_end_z = sample_df['prot_epz'].values   # конец - детектирование протона
    prot_end_x = sample_df['prot_epx'].values
    
    # Координаты для пи-минусов (зеленые) - тоже от точки распада лямбды
    pimin_start_z = sample_df['lam_epz'].values  # начало - точка распада лямбды
    pimin_start_x = sample_df['lam_epx'].values
    pimin_end_z = sample_df['pimin_epz'].values  # конец - детектирование пиона
    pimin_end_x = sample_df['pimin_epx'].values
    
    # Создаем сегменты для каждой частицы
    lam_segments = np.array([[[lam_start_z[i], lam_start_x[i]], [lam_end_z[i], lam_end_x[i]]] for i in range(len(sample_df))])
    prot_segments = np.array([[[prot_start_z[i], prot_start_x[i]], [prot_end_z[i], prot_end_x[i]]] for i in range(len(sample_df))])
    pimin_segments = np.array([[[pimin_start_z[i], pimin_start_x[i]], [pimin_end_z[i], pimin_end_x[i]]] for i in range(len(sample_df))])
    
    # Создаем график
    fig, ax = create_plot_with_background(bck_image="eic_center_forward_bw.png")
    
    # Добавляем линии траекторий
    lam_lines = LineCollection(lam_segments, color='red', alpha=0.7, linewidths=2, label='Λ⁰ trajectory')
    prot_lines = LineCollection(prot_segments, color='blue', alpha=0.7, linewidths=2, label='Proton trajectory')  
    pimin_lines = LineCollection(pimin_segments, color='lime', alpha=0.7, linewidths=2, label='π⁻ trajectory')
    
    ax.add_collection(lam_lines)
    ax.add_collection(prot_lines)
    ax.add_collection(pimin_lines)
    
    # Добавляем точки для наглядности
    ax.scatter(lam_start_z, lam_start_x, color='darkred', s=50, marker='o', label='Λ⁰ birth', alpha=1)
    ax.scatter(lam_end_z, lam_end_x, color='red', s=50, marker='s', label='Λ⁰ decay', alpha=1)
    ax.scatter(prot_end_z, prot_end_x, color='blue', s=50, marker='^', label='Proton detection', alpha=1)
    ax.scatter(pimin_end_z, pimin_end_x, color='lime', s=50, marker='v', label='π⁻ detection', alpha=1)
    
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title("Λ⁰ → p + π⁻ decay trajectories")
    
    plt.tight_layout()
    plt.show()



def plot_neutron_pizero_decay_trajectories():
        # взять только 20 евентов # (столько, сколько хорошо видно разницу) построить так: красная паочка как летят лямбды (точка старт - ф ниш),
    #  синяя палочка - как летят протоны, зеленая - пионы

    sample_df = n_pi_zero_decays.iloc[34:200].head(1)
    # Координаты для лямбд (красные)
    lam_start_z = sample_df['lam_vz'].values  # начало - рождение
    lam_start_x = sample_df['lam_vx'].values
    lam_end_z = sample_df['lam_epz'].values   # конец - распад  
    lam_end_x = sample_df['lam_epx'].values

    # Координаты для протонов (синие) - летят от точки распада лямбды
    prot_start_z = sample_df['lam_epz'].values  # начало - точка распада лямбды
    prot_start_x = sample_df['lam_epx'].values  
    prot_end_z = sample_df['neut_epz'].values   # конец - детектирование протона
    prot_end_x = sample_df['neut_epx'].values

    # Координаты для пи-минусов (зеленые) - тоже от точки распада лямбды
    pimin_start_z = sample_df['lam_epz'].values  # начало - точка распада лямбды
    pimin_start_x = sample_df['lam_epx'].values
    pimin_end_z = sample_df['pizero_epz'].values  # конец - детектирование пиона
    pimin_end_x = sample_df['pizero_epx'].values

    # Создаем сегменты для каждой частицы
    lam_segments = np.array([[[lam_start_z[i], lam_start_x[i]], [lam_end_z[i], lam_end_x[i]]] for i in range(len(sample_df))])
    prot_segments = np.array([[[prot_start_z[i], prot_start_x[i]], [prot_end_z[i], prot_end_x[i]]] for i in range(len(sample_df))])
    pimin_segments = np.array([[[pimin_start_z[i], pimin_start_x[i]], [pimin_end_z[i], pimin_end_x[i]]] for i in range(len(sample_df))])

    # Создаем график
    fig, ax = create_plot_with_background(bck_image="eic_center_forward_bw.png")

    # Добавляем линии траекторий
    lam_lines = LineCollection(lam_segments, color='red', alpha=0.7, linewidths=2, label='Λ⁰ trajectory')
    prot_lines = LineCollection(prot_segments, color='blue', alpha=0.7, linewidths=2, label='Proton trajectory')  
    pimin_lines = LineCollection(pimin_segments, color='lime', alpha=0.7, linewidths=2, label='π⁻ trajectory')

    ax.add_collection(lam_lines)
    ax.add_collection(prot_lines)
    ax.add_collection(pimin_lines)

    # Добавляем точки для наглядности
    ax.scatter(lam_start_z, lam_start_x, color='darkred', s=50, marker='o', label='Λ⁰ birth', alpha=1)
    ax.scatter(lam_end_z, lam_end_x, color='red', s=50, marker='s', label='Λ⁰ decay', alpha=1)
    ax.scatter(prot_end_z, prot_end_x, color='blue', s=50, marker='^', label='Proton detection', alpha=1)
    ax.scatter(pimin_end_z, pimin_end_x, color='lime', s=50, marker='v', label='π⁻ detection', alpha=1)

    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title("Λ⁰ → p + π⁻ decay trajectories")

    plt.tight_layout()
    plt.show()


def plot_undecayed_primary_lambdas():
    undecayed_primary_percentage = calculate_percentage(df, df['lam_is_first'] == 0)

    print(f"Undecayed primary lambdas: {undecayed_primary_percentage:.1f}%")
    if undecayed_primary_percentage == 0:
        plot_point(df_not_dec['lam_epz'], df_not_dec['lam_epx'])

   
def plot_primary_vs_secondary_decay_points():
    """График где распадаются праймари и не праймари лямбды (точки)"""
    
    fig, ax = create_plot_with_background()
    
    ax.scatter(df_pram['lam_epz'], df_pram['lam_epx'], 
               color='red', marker="o", s=20, alpha=0.5, 
               label=f'Primary Λ⁰')
    
    ax.scatter(df_not_pram['lam_epz'], df_not_pram['lam_epx'], 
               color='blue', marker="s", s=20, alpha=0.5,
               label=f'Secondary Λ⁰')
    
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_title("Primary vs Secondary Λ⁰ decay points")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.show()


def plot_primary_vs_secondary_birth_points():
    """График где рождаются праймари и не праймари лямбды (точки)"""
    
    fig, ax = create_plot_with_background()
    
    ax.scatter(df_not_pram['lam_vz'], df_not_pram['lam_vx'], 
                   color='blue', marker="s", s=20, alpha=0.5,
                   label=f'Secondary Λ⁰')

    ax.scatter(df_pram['lam_vz'], df_pram['lam_vx'], 
               color='red', marker="o", s=20, alpha=0.5, 
               label=f'Primary Λ⁰')
    
    
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_title("Primary vs Secondary Λ⁰ birth points")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.show()

def main():

    """
    Сколько лямбд которые распались на протон или нейтрон или не распались вообще. (проценты)

    Есть ли лямбды, которые праймари и не распались вообще? Если они есть, показать на графике, где они распадаются.

    Сколько лямбд, которые не праймари? Сколько эвентов со вторичными лямбдами? (все распады)

    График где распадаются праймари и не праймари лямбды (точки)

    График где рождаются праймари и не праймаи лямбды (точки)

    1д гистограмма зет точки распада лямбды (только праймари)

    2Д гистограмма на фоне детектора, где распадаются лямбды.
    --
    """
    calculate_decayed() 
    plot_undecayed_primary_lambdas()

    # Доделать
    analyze_primary_vs_secondary_lambdas()

    plot_primary_vs_secondary_decay_points()

    plot_primary_vs_secondary_birth_points()

    plot_primary_lambda_decay_z_distribution()

    plt_hist2d(df_pram['lam_epz'], df_pram['lam_epx'], 
           title="Λ⁰ decay points distribution")


    """
    2Д гистограмма на фоне детектора, где распадаются лямбды (for proton+pi-)

    взять только 20 евентов (столько, сколько хорошо видно разницу) построить так:
      красная паочка как летят лямбды (точка старт - ф ниш), синяя палочка - как летят протоны, зеленая - пионы

    2Д гистограмма пролета протона

    2Д гистограмма пролета пиона

    2Д гистограмма точки конца протона

    2Д гистограмма точки контца пиона
    """

    plt_hist2d(p_pi_minus_decays['lam_epz'], p_pi_minus_decays['lam_epx'], 
           title="Λ⁰ → p + π⁻ decay points distribution")
    
    decay_trajectories()

    plt_hist2d(p_pi_minus_decays['pimin_epz'], p_pi_minus_decays['pimin_epx'], 
           title="π⁻ end points distribution")


    """

    2Д гистограмма на фоне детектора, где распадаются лямбды (for neutron + pi0)

    взять только 20 евентов (столько, сколько хорошо видно разницу) построить так: красная паочка как летят лямбды (точка старт - ф ниш), синяя палочка - как летят протоны, зеленая - пионы

    2Д гистограмма пролета нейтрона (-)

    2Д гистограмма пролета пиона (-)

    конец нейтрона (!)

    2Д гистограмма точки конца/разавла пи0 (!)

    2Д гистограмма точки конца gamma gamma (!)
    """

    plt_hist2d(n_pi_zero_decays['lam_epz'], n_pi_zero_decays['lam_epx'], 
           title="n + π⁰ decay points distribution")
    
    plot_neutron_pizero_decay_trajectories()


    plt_hist2d(n_pi_zero_decays['pizero_epz'], n_pi_zero_decays['pizero_epx'],
           title="π⁰ decay points distribution")
    
    # 2Д гистограмма точки конца gamma gamma
    plt_hist2d(n_pi_zero_decays['gamone_epz'], n_pi_zero_decays['gamone_epx'],bins=150,
           title="Gamma end points distribution")


    


if __name__ == "__main__":
    main()
