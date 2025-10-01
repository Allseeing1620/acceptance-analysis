
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
].copy()

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


def plt_hist2d(x_axis, y_axis, bins=50, bin_size=None, 
               xlabel="z [mm]", ylabel="y [mm]", title="Title is missing"):

    fig, ax = create_plot_with_background()
    
    original_xlim = ax.get_xlim()
    original_ylim = ax.get_ylim()
    
    if bin_size is not None:
        x_range = original_xlim[1] - original_xlim[0]
        y_range = original_ylim[1] - original_ylim[0]
        
        x_bins = max(1, int(x_range / bin_size))
        y_bins = max(1, int(y_range / bin_size))
        bins = [x_bins, y_bins]
        
        print(f"Размер бина: {bin_size}×{bin_size} мм")
        print(f"Количество бинов: {x_bins}×{y_bins}")
    
    h = ax.hist2d(x_axis, y_axis, bins=bins,
                   cmap='viridis', norm=LogNorm())

    ax.set_xlim(original_xlim)
    ax.set_ylim(original_ylim)
    
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

    sample_df = n_pi_zero_decays.iloc[34:200].head(1)
    
    # Координаты для лямбд (красные)
    lam_start_z = sample_df['lam_vz'].values
    lam_start_x = sample_df['lam_vx'].values
    lam_end_z = sample_df['lam_epz'].values
    lam_end_x = sample_df['lam_epx'].values

    # Координаты для нейтронов (синие)
    neut_start_z = sample_df['lam_epz'].values
    neut_start_x = sample_df['lam_epx'].values  
    neut_end_z = sample_df['neut_epz'].values
    neut_end_x = sample_df['neut_epx'].values

    # Координаты для пи-нолей (зеленые)
    pizero_start_z = sample_df['lam_epz'].values
    pizero_start_x = sample_df['lam_epx'].values
    pizero_end_z = sample_df['pizero_epz'].values
    pizero_end_x = sample_df['pizero_epx'].values

    # Создаем сегменты
    lam_segments = np.array([[[lam_start_z[i], lam_start_x[i]], [lam_end_z[i], lam_end_x[i]]] for i in range(len(sample_df))])
    neut_segments = np.array([[[neut_start_z[i], neut_start_x[i]], [neut_end_z[i], neut_end_x[i]]] for i in range(len(sample_df))])
    pizero_segments = np.array([[[pizero_start_z[i], pizero_start_x[i]], [pizero_end_z[i], pizero_end_x[i]]] for i in range(len(sample_df))])

    # Создаем график
    fig, ax = create_plot_with_background(bck_image="eic_center_forward_bw.png")

    # Добавляем линии траекторий
    lam_lines = LineCollection(lam_segments, color='red', alpha=0.7, linewidths=2, label='Λ⁰ trajectory')
    neut_lines = LineCollection(neut_segments, color='blue', alpha=0.7, linewidths=2, label='Neutron trajectory')  
    pizero_lines = LineCollection(pizero_segments, color='lime', alpha=0.7, linewidths=2, label='π⁰ trajectory')

    ax.add_collection(lam_lines)
    ax.add_collection(neut_lines)
    ax.add_collection(pizero_lines)

    # Добавляем точки
    ax.scatter(lam_start_z, lam_start_x, color='darkred', s=50, marker='o', label='Λ⁰ birth', alpha=1)
    ax.scatter(lam_end_z, lam_end_x, color='red', s=50, marker='s', label='Λ⁰ decay', alpha=1)
    ax.scatter(neut_end_z, neut_end_x, color='blue', s=50, marker='^', label='Neutron detection', alpha=1)
    ax.scatter(pizero_end_z, pizero_end_x, color='lime', s=50, marker='v', label='π⁰ detection', alpha=1)

    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title("Neutron + π⁰ decay trajectories")

    plt.tight_layout()
    plt.show()


def plot_undecayed_primary_lambdas():
    undecayed_primary = df_pram[(df_pram['prot_id'].isna() & df_pram['neut_id'].isna())]
    undecayed_primary_percentage = calculate_percentage(df, undecayed_primary.index.isin(df.index))

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

import os
from datetime import datetime

def save_all_plots(output_dir="analysis_results"):
    """
    Сохраняет все графики в указанную папку
    """
    import os
    from datetime import datetime
    
    # Создаем папку если не существует
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = os.path.join(output_dir, f"analysis_{timestamp}")
    os.makedirs(results_dir, exist_ok=True)
    
    print(f"Saving plots to: {results_dir}")
    
    def save_plot(func, filename, *args, **kwargs):
        try:
            plt.switch_backend('Agg')
            func(*args, **kwargs)
            plt.savefig(os.path.join(results_dir, filename), dpi=150, bbox_inches='tight')
            plt.close()
            print(f"✓ Saved: {filename}")
        except Exception as e:
            print(f"✗ Error saving {filename}: {e}")
        finally:
            plt.switch_backend('TkAgg')  # возвращаем нормальный backend
    
    # Сохраняем ВСЕ графики которые есть в main()
    save_plot(calculate_decayed, "01_decay_statistics.png")
    save_plot(plot_undecayed_primary_lambdas, "02_undecayed_primary_lambdas.png")
    save_plot(plot_primary_vs_secondary_decay_points, "03_primary_vs_secondary_decay_points.png")
    save_plot(plot_primary_vs_secondary_birth_points, "04_primary_vs_secondary_birth_points.png")
    save_plot(plot_primary_lambda_decay_z_distribution, "05_primary_lambda_decay_z_distribution.png")
    
    # 2D гистограммы через lambda
    save_plot(lambda: plt_hist2d(df_pram['lam_epz'], df_pram['lam_epx'], 
                                title="Λ⁰ decay points distribution"), 
             "06_lambda_decay_points.png")
    
    save_plot(lambda: plt_hist2d(p_pi_minus_decays['lam_epz'], p_pi_minus_decays['lam_epx'], 
                                title="Λ⁰ → p + π⁻ decay points"), 
             "07_proton_pion_decay_points.png")
    
    save_plot(decay_trajectories, "08_proton_pion_trajectories.png")
    
    save_plot(lambda: plt_hist2d(p_pi_minus_decays['pimin_epz'], p_pi_minus_decays['pimin_epx'], 
                                title="π⁻ end points"), 
             "09_pion_end_points.png")
    
    save_plot(lambda: plt_hist2d(n_pi_zero_decays['lam_epz'], n_pi_zero_decays['lam_epx'], 
                                title="n + π⁰ decay points"), 
             "10_neutron_pizero_decay_points.png")
    
    save_plot(plot_neutron_pizero_decay_trajectories, "11_neutron_pizero_trajectories.png")
    
    save_plot(lambda: plt_hist2d(n_pi_zero_decays['pizero_epz'], n_pi_zero_decays['pizero_epx'],
                                title="π⁰ decay points"), 
             "12_pizero_decay_points.png")
    
    save_plot(lambda: plt_hist2d(n_pi_zero_decays['gamone_epz'], n_pi_zero_decays['gamone_epx'], bins=150,
                                title="Gamma end points"), 
             "13_gamma_end_points.png")
    
    print(f"\n All plots saved to: {results_dir}")
    return results_dir

def plot_particle_trajectory_histogram(particle_type, dataframe, 
                                     start_x_col, start_z_col,
                                     end_x_col, end_z_col,
                                     grid_x_step=100, grid_z_step=100,
                                     min_trajectories=10, cmap='viridis',
                                     particle_name=None):
    """
    Универсальная функция для построения 2D гистограммы траекторий частиц
    
    Parameters:
    -----------
    particle_type : str
        Тип частицы ('proton', 'pion', 'kaon', etc.) - для названий
    dataframe : DataFrame
        DataFrame с данными о частицах
    start_x_col, start_z_col : str
        Названия колонок с координатами начала траектории
    end_x_col, end_z_col : str
        Названия колонок с координатами конца траектории  
    grid_x_step, grid_z_step : float
        Размер ячейки сетки в мм
    min_trajectories : int
        Минимальное количество траекторий для отображения
    cmap : str
        Цветовая карта
    particle_name : str, optional
        Отображаемое имя частицы (если None, берется из particle_type)
    """
    
    if particle_name is None:
        particle_name = particle_type
    
    fig, ax = create_plot_with_background()
    GRID_Z_MIN, GRID_Z_MAX = ax.get_xlim()
    GRID_X_MIN, GRID_X_MAX = ax.get_ylim()
    
    original_xticks = ax.get_xticks()
    original_yticks = ax.get_yticks()
    
    plt.close(fig)
    
    GRID_X_BINS = int((GRID_X_MAX - GRID_X_MIN) / grid_x_step)
    GRID_Z_BINS = int((GRID_Z_MAX - GRID_Z_MIN) / grid_z_step)
    
    
    
    def get_grid_index(x, z):
        x_idx = int((x - GRID_X_MIN) / grid_x_step)
        z_idx = int((z - GRID_Z_MIN) / grid_z_step)
        x_idx = max(0, min(GRID_X_BINS - 1, x_idx))
        z_idx = max(0, min(GRID_Z_BINS - 1, z_idx))
        return (x_idx, z_idx)
    
    def find_intersected_cells_bresenham(start, end):
        start_x, start_z = start
        end_x, end_z = end
        
        x0_idx, z0_idx = get_grid_index(start_x, start_z)
        x1_idx, z1_idx = get_grid_index(end_x, end_z)
        
        cells = []
        
        dx = abs(x1_idx - x0_idx)
        dy = abs(z1_idx - z0_idx)
        
        x, z = x0_idx, z0_idx
        
        x_step = 1 if x1_idx > x0_idx else -1
        z_step = 1 if z1_idx > z0_idx else -1
        
        if dx > dy:
            err = dx / 2.0
            while x != x1_idx:
                cells.append((x, z))
                err -= dy
                if err < 0:
                    z += z_step
                    err += dx
                x += x_step
        else:
            err = dy / 2.0
            while z != z1_idx:
                cells.append((x, z))
                err -= dx
                if err < 0:
                    x += x_step
                    err += dy
                z += z_step
        
        cells.append((x1_idx, z1_idx))
        return list(set(cells))
    
    
    grid_hist = np.zeros((GRID_Z_BINS, GRID_X_BINS))
    
    
    for i in range(len(dataframe)):
        particle = dataframe.iloc[i]
        
        
        start_point = (particle[start_x_col], particle[start_z_col])
        end_point = (particle[end_x_col], particle[end_z_col])
        
    
        intersected_cells = find_intersected_cells_bresenham(start_point, end_point)
        
       
        for x_idx, z_idx in intersected_cells:
            grid_hist[z_idx, x_idx] += 1
    
    
    fig, ax = create_plot_with_background()
    ax.grid(False)
    
    
    x_edges = np.linspace(GRID_X_MIN, GRID_X_MAX, GRID_X_BINS + 1)
    z_edges = np.linspace(GRID_Z_MIN, GRID_Z_MAX, GRID_Z_BINS + 1)
    
    
    grid_hist_masked = np.where(grid_hist >= min_trajectories, grid_hist, np.nan)
    
    
    im = ax.pcolormesh(z_edges, x_edges, grid_hist_masked.T, 
                      cmap=cmap, alpha=0.8, shading='auto')
    
    
    cbar = plt.colorbar(im, ax=ax, 
                       label=f'Number of {particle_name} trajectories (min {min_trajectories})', 
                       shrink=0.6, aspect=20, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x [mm]")
    ax.set_title(f"{particle_name.capitalize()} Trajectory Histogram\n(Density of {particle_name} paths through detector)")
    ax.set_aspect("equal", adjustable="box")
    
    
    ax.set_xticks(original_xticks)
    ax.set_yticks(original_yticks)
    
    ax.set_xticks(np.arange(GRID_Z_MIN, GRID_Z_MAX + grid_z_step, grid_z_step), minor=True)
    ax.set_yticks(np.arange(GRID_X_MIN, GRID_X_MAX + grid_x_step, grid_x_step), minor=True)
    ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, which='minor')
    ax.grid(True, alpha=0.4, linestyle='-', linewidth=1, which='major')
    
    plt.tight_layout()
    plt.show()
    
    return grid_hist
def main():

    """
    Сколько лямбд которые распались на протон или нейтрон или не распались вообще. (проценты) 

    Есть ли лямбды, которые праймари и не распались вообще? Если они есть, показать на графике, где они распадаются. 

    Сколько лямбд, которые не праймари? Сколько эвентов со вторичными лямбдами? (все распады) ДОБАВИТЬ ГИСТОГРАММЫ СКОЛЬКО ЛЯМБД В СОБЫТИИ посмотреть. 
    добавить может 0 тк есть лямбы без лямбд

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
    
    plot_particle_trajectory_histogram(
    particle_type='proton',
    dataframe=p_pi_minus_decays,
    start_x_col='prot_vx', start_z_col='prot_vz',
    end_x_col='prot_epx', end_z_col='prot_epz',
    grid_x_step=100, grid_z_step=100,
    min_trajectories=10,
    cmap='viridis'
    )

    plot_particle_trajectory_histogram(
    particle_type='pion',
    dataframe=p_pi_minus_decays, 
    start_x_col='pimin_vx', start_z_col='pimin_vz',
    end_x_col='pimin_epx', end_z_col='pimin_epz',
    grid_x_step=100, grid_z_step=100,
    min_trajectories=100,
    cmap='viridis'
    )

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
    
    

    plot_particle_trajectory_histogram(
    particle_type='neut',
    dataframe=n_pi_zero_decays,
    start_x_col='neut_vx', start_z_col='neut_vz',
    end_x_col='neut_epx', end_z_col='neut_epz',
    grid_x_step=100, grid_z_step=100,
    min_trajectories=10,
    cmap='viridis'
    )

    plot_particle_trajectory_histogram(
    particle_type='pizero',
    dataframe=n_pi_zero_decays, 
    start_x_col='pizero_vx', start_z_col='pizero_vz',
    end_x_col='pizero_epx', end_z_col='pizero_epz',
    grid_x_step=300, grid_z_step=300,
    min_trajectories=100,
    cmap='viridis'
    )




    #save_all_plots("C:/Users/User/Desktop/test/my_analysis_results")


if __name__ == "__main__":
    main()
