import os
import pandas as pd

def modify_hydrotable_bathymetric_adjustments(main_directory, huc8, df_bfull_geom):
    """
    Modify the hydrotable by applying bathymetric adjustments.

    Parameters:
        main_directory (str): The base directory containing flood data.
        huc8 (str): Hydrologic Unit Code.
        df_bfull_geom (pd.DataFrame): DataFrame containing bathymetric information.

    Returns:
        pd.DataFrame: The modified hydrotable DataFrame.
    """
    branches_path = main_directory + f"/flood_{huc8}/{huc8}/branches"
    folders = [f for f in os.listdir(branches_path) if os.path.isdir(os.path.join(branches_path, f))]
    
    hydrotable = pd.DataFrame()

    for branch in folders:
        src_filename = branches_path+f'/{branch}/src_full_crosswalked_{branch}.csv'
                                    
        src_df = pd.read_csv(src_filename)
        src_df = src_df.rename(columns={'Discharge (m3s-1)': 'Discharge (m3s-1)_Org'})
        src_df = pd.merge(src_df, df_bfull_geom, how='left', on='feature_id')
        
        src_df = src_df.rename(columns={
            'bnk_width': 'top_width',
            'bnk_xsce_A': 'missing_xs_area_m2',
            'bnk_wet_p': 'missing_wet_perimeter_m',
            'mf_depth': 'shift_datume_m'
        })
        
        src_df['Bathymetry_source'] = "Zarrabi_etal"
        src_df["BARC"] = "TRUE"
        
        # Fill missing values
        src_df[['missing_xs_area_m2', 'missing_wet_perimeter_m']] = src_df[
            ['missing_xs_area_m2', 'missing_wet_perimeter_m']
        ].fillna(0.0)

        # Compute adjusted Volume and Bed Area
        src_df.loc[src_df["Bathymetry_source"] == "Zarrabi_etal", "Volume (m3)"] = (
            src_df['Volume (m3)'] + (src_df['missing_xs_area_m2'] * (src_df['LENGTHKM'] * 1000))
        )
        
        src_df.loc[src_df["Bathymetry_source"] == "Zarrabi_etal", "BedArea (m2)"] = (
            src_df['BedArea (m2)'] + (src_df['missing_wet_perimeter_m'] * (src_df['LENGTHKM'] * 1000))
        )

        # Recalculate hydraulic properties
        src_df.loc[src_df["Bathymetry_source"] == "Zarrabi_etal", "WettedPerimeter (m)"] = (
            src_df['WettedPerimeter (m)'] + src_df['missing_wet_perimeter_m']
        )
        
        src_df.loc[src_df["Bathymetry_source"] == "Zarrabi_etal", "WetArea (m2)"] = (
            src_df['WetArea (m2)'] + src_df['missing_xs_area_m2']
        )

        # Compute hydraulic radius
        src_df["HydraulicRadius (m)"] = src_df["WetArea (m2)"] / src_df["WettedPerimeter (m)"]
        src_df["HydraulicRadius (m)"] = src_df["HydraulicRadius (m)"].fillna(0)

        # Compute new discharge values
        src_df["Discharge (m3s-1)_BARC"] = (
            src_df["WetArea (m2)"] * 
            (src_df["HydraulicRadius (m)"] ** (2.0 / 3)) * 
            (src_df["SLOPE"] ** 0.5) / 
            src_df["ManningN"]
        )
        
        # Save updated src file
        src_df.to_csv(src_filename, index=False)

        # Prepare trimmed DataFrame for merging
        df_src_trim = src_df[["HydroID", "Stage", "Bathymetry_source", "Discharge (m3s-1)_BARC", "shift_datume_m"]]
        df_src_trim = df_src_trim.rename(columns={"Stage": "stage"})
        df_src_trim["discharge_cms"] = df_src_trim["Discharge (m3s-1)_BARC"]

        # Read and update hydroTable
        htable_filename = branches_path+f'/{branch}/hydroTable_{branch}.csv'
        df_htable = pd.read_csv(htable_filename, dtype={"HUC": str, "last_updated": object, "submitter": object, "obs_source": object})
        df_htable["discharge_cms_org"] = df_htable["discharge_cms"]

        # Drop previous discharge-related columns
        df_htable = df_htable.drop(
            [
                "subdiv_applied",
                "discharge_cms",
                "overbank_n",
                "channel_n",
                "subdiv_discharge_cms",
                "Bathymetry_source",
            ],
            axis=1,
            errors="ignore",
        )

        # Merge updated data into hydroTable
        df_htable = df_htable.merge(
            df_src_trim, how="left", left_on=["HydroID", "stage"], right_on=["HydroID", "stage"]
        )

        htable= pd.DataFrame()
        
        for i in df_htable['HydroID'].unique():
            df_plt = df_htable[df_htable['HydroID'] == i]

            # Adjust stage values
            df_plt['stage'] = df_plt['stage'] + df_plt['shift_datume_m']

            # Add a zero-value row
            first_row = df_plt.iloc[0]
            first_row[['stage', 'discharge_cms']] = 0 

            # Insert at the beginning
            df_plt = pd.concat([pd.DataFrame([first_row]), df_plt], ignore_index=True)

            # Append to the final DataFrame
            htable = pd.concat([htable, df_plt], ignore_index=True)

        

        # Save updated hydroTable file
        htable.to_csv(htable_filename, index=False)

        # Append branch ID for tracking
        htable["branch_id"] = branch
        hydrotable = pd.concat([hydrotable, htable])

    # Move 'branch_id' to the second column
    cols = list(hydrotable.columns)
    cols.insert(1, cols.pop(cols.index("branch_id")))
    hydrotable = hydrotable[cols]

    # Save final hydrotable file
    output_filename = main_directory + f"/flood_{huc8}/{huc8}/hydrotable.csv"
    hydrotable.to_csv(output_filename, index=False)

    return hydrotable

