import pandas as pd
import numpy as np


def get_inputs(period, dt):
    """
    Returns input data as in a DataFrame, where:
    index - DateTime
    columns - input data

    E.g. if period = 8 h and dt = 60 minutes, it will return
    a DataFrame with 9 rows, starting from t = Jan 1, 00:00:00 and ending at
    t = Jan 1, 08:00:00.

    :param period: int, how much data to read [h]
    :param period: int, sampling rate [minutes]
    :return: DataFrame
    """

    # Time steps as in the CSV files
    t_start = 0   # [s]
    t_final = period * 3600  # [s]
    csv_dt = 60.  # [s]
    n_steps = int((t_final - t_start) / csv_dt) + 1

    # Read CSV files
    Tout = pd.read_csv('resources/ou44/data/df_tamb.csv')['weather'].iloc[:n_steps].values
    Hglo = pd.read_csv('resources/ou44/data/df_solrad.csv')['weather'].iloc[:n_steps].values
    qhvac = np.full(n_steps, 0.)
    nocc = pd.read_csv('resources/ou44/data/df_occ.csv')['22-508-1'].iloc[:n_steps].values

    # Put into DataFrame, resample and extract again
    index = pd.to_datetime(
        pd.read_csv('resources/ou44/data/df_tamb.csv')['datetime'][:n_steps]
    )
    df = pd.DataFrame(index=index)
    df['Tout'] = Tout
    df['Hglo'] = Hglo
    df['qhvac'] = qhvac
    df['nocc'] = nocc
    df = df.resample('{}min'.format(dt)).mean()
    Tout = df['Tout'].values
    Hglo = df['Hglo'].values
    qhvac = df['qhvac'].values
    nocc = df['nocc'].values

    # u, t, u_names
    # t = np.arange(t_start, t_final, dt)
    # u = np.vstack((Tout, Hglo, qhvac, nocc)).T
    # u_names = ['Tout', 'Hglo', 'qhvac', 'nocc']

    return df
