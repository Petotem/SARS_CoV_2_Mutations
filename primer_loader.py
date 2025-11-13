import pandas as pd
import numpy as np


class Primer:
    """A class to represent primer information."""

    def __init__(self):
        self.id = ""
        self.primers = []
        self.lengths = []
        self.ids = []
        self.start_position = 0
        self.nsp = ""
        self.group = ""
        self.subgroup = ""
        self.type = ""
        self.temperature = 0.0

        self.align_score = 0
        self.exact_match = 0
        self.partial_match = 0
        self.total_sequence = 0


class PrimerOperation:
    """Handles loading and grouping primer data from Excel files."""

    def __init__(self):
        self.primer_dict = {}

    def load_primers_from_file(self, filename: str):
        """
        Load primer information from an Excel file and structure it by group and subgroup.

        Parameters
        ----------
        filename : str
            Path to the Excel file containing primer data.

        Returns
        -------
        tuple
            primer_dict: dict
                A dictionary mapping each group to a list of Primer objects.
            all_primers: list
                A list of all primer IDs.
        """
        df = pd.read_excel(filename)

        primer_dict = {}
        groups = df.groupby("Group")

        for group_name, group_df in groups:
            primer_info_list = []

            for subgroup_name, subgroup_df in group_df.groupby("SubGroup"):
                primers_array = subgroup_df.to_numpy()

                primer_obj = Primer()
                primer_obj.group = primers_array[0, df.columns.get_loc("Group")]
                primer_obj.subgroup = primers_array[0, df.columns.get_loc("SubGroup")]
                primer_obj.nsp = primers_array[0, df.columns.get_loc("NSP")] if "NSP" in df.columns else ""
                primer_obj.start_position = primers_array[0, df.columns.get_loc("Position")] if "Position" in df.columns else 0
                primer_obj.type = primers_array[0, df.columns.get_loc("Type")] if "Type" in df.columns else ""

                # Clean and convert temperature
                temp_raw = str(primers_array[0, df.columns.get_loc("Temperature")]).replace('\t', '').replace('ºC', '').replace(' ', '')
                try:
                    primer_obj.temperature = float(temp_raw)
                except ValueError:
                    primer_obj.temperature = np.nan

                # Extract sequences and IDs
                primer_obj.primers = subgroup_df["PrimerSequence"].tolist() if "PrimerSequence" in df.columns else []
                primer_obj.lengths = [len(p) for p in primer_obj.primers]
                primer_obj.ids = subgroup_df.iloc[:, 0].tolist()

                primer_info_list.append(primer_obj)

            primer_dict[group_name] = primer_info_list

        all_primers = df.iloc[:, 0].tolist()

        return primer_dict, all_primers
