import pandas as pd
from pathlib import Path

class Sequence:
    """Represents a single genome sequence record."""
    
    VALID_LENGTH = 29903
    VALID_RANGE = 500
    VALID_BASES = {'A', 'T', 'G', 'C'}
    SEPARATOR = "\t"

    def __init__(self):
        self.id = ""
        self.header = ""
        self.epi_id = "NA"
        self.date = "NA"
        self.sequence = ""
        self.country = "NA"
        self.continent = "NA"
        self.length = 0
        self.age = "NA"
        self.sex = "NA"
        self.host = "NA"
        self.submitted_date = "NA"
        self.virus_type = "NA"
        self.lineage = "NA"
        self.error = False

    def is_length_valid(self) -> bool:
        return (self.VALID_LENGTH - self.VALID_RANGE) <= self.length <= (self.VALID_LENGTH + self.VALID_RANGE)

    def is_human(self) -> bool:
        return self.host.lower() == "human"

    def contains_invalid_bases(self) -> bool:
        return any(base not in self.VALID_BASES for base in self.sequence)

    @classmethod
    def from_line(cls, line: str):
        seq = cls()
        try:
            parts = line.strip().split(cls.SEPARATOR)
            seq.country = parts[0]
            seq.sequence = parts[1]
            seq.length = int(parts[2])
            seq.id = parts[3]
            seq.epi_id = parts[4]
            seq.date = parts[5]
            seq.continent = parts[6]
            seq.sex = parts[7]
            seq.age = parts[8]
            seq.host = parts[9]
        except Exception:
            seq.error = True
        return seq


class CovidFastaProcessor:
    """Processes a FASTA file and merges with metadata for COVID-19 primer analysis."""

    def __init__(self, fasta_path: str, metadata_path: str = None, lineage_path: str = None,
                 output_dir: str = "result", split_limit: int = 5_000_000):
        self.fasta_path = Path(fasta_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.split_limit = split_limit

        self.metadata_df = None
        self.lineage_df = None

        if metadata_path:
            self.metadata_df = pd.read_csv(metadata_path, sep="\t", index_col="Accession ID")
        if lineage_path:
            self.lineage_df = pd.read_excel(lineage_path)

    def _is_valid_sequence(self, seq: Sequence) -> bool:
        return not seq.error and seq.is_length_valid()

    def _split_header(self, header_line: str) -> Sequence:
        seq = Sequence()
        try:
            parts = header_line.strip().split("|")
            seq.header = parts[0].replace(">", "")
            seq.id = parts[1]
            seq.epi_id = parts[1]
            seq.date = parts[2] if len(parts) > 2 else "NA"
        except Exception:
            seq.error = True
        return seq

    def _read_sequence(self, file_handle) -> tuple[str, str]:
        """Read one sequence block from FASTA."""
        sequence_lines = []
        while True:
            line = file_handle.readline()
            if not line or line.startswith(">"):
                return "".join(sequence_lines).replace("\n", ""), line
            sequence_lines.append(line.strip())

    def _find_lineage(self, lineage: str) -> str:
        if self.lineage_df is None:
            return "Unknown"
        for _, row in self.lineage_df.iterrows():
            if str(lineage).lower().find(str(row[1]).lower()) >= 0:
                return row[0]
        return "Other"

    def _attach_metadata(self, seq: Sequence) -> Sequence:
        if self.metadata_df is None:
            return seq
        try:
            row = self.metadata_df.loc[seq.epi_id]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            seq.sex = row.get("Gender", "NA")
            location = str(row.get("Location", "NA")).split("/")
            seq.country = location[1].strip() if len(location) > 1 else "NA"
            seq.date = row.get("Collection date", "NA")
            seq.lineage = self._find_lineage(row.get("Lineage", ""))
        except KeyError:
            pass
        return seq

    def _save_sequence(self, seq: Sequence, file_handle):
        file_handle.write(
            f"{seq.country},{seq.sequence},{seq.length},{seq.id},{seq.epi_id},{seq.date},"
            f"{seq.submitted_date},{seq.continent},{seq.sex},{seq.age},{seq.host},"
            f"{seq.virus_type},{seq.lineage},{seq.header}\n"
        )

    def process(self, output_prefix: str = "split", stats_file: str = "statistic"):
        total, valid, count, file_idx = 0, 0, 0, 0
        output_file = self.output_dir / f"{output_prefix}_{file_idx}.csv"

        with open(self.fasta_path, "r") as reader, open(output_file, "w") as writer:
            writer.write("Country,Sequence,Length,Id,EpiID,Date,SubmitDate,Continent,Sex,Age,Host,VirusType,Lineage,Header\n")

            line = reader.readline()
            while line:
                if not line.startswith(">"):
                    line = reader.readline()
                    continue

                seq_info = self._split_header(line)
                sequence, next_line = self._read_sequence(reader)
                seq_info.sequence = sequence
                seq_info.length = len(sequence)
                total += 1
                seq_info = self._attach_metadata(seq_info)

                if self._is_valid_sequence(seq_info):
                    valid += 1
                    self._save_sequence(seq_info, writer)
                    count += 1
                    if count >= self.split_limit:
                        file_idx += 1
                        count = 0
                        writer.close()
                        writer = open(self.output_dir / f"{output_prefix}_{file_idx}.csv", "w")
                        writer.write("Country,Sequence,Length,Id,EpiID,Date,SubmitDate,Continent,Sex,Age,Host,VirusType,Lineage,Header\n")

                line = next_line

        stats_path = self.output_dir / f"{stats_file}.csv"
        with open(stats_path, "w") as stats:
            stats.write(f"Total: {total}\nValid: {valid}\n")


if __name__ == "__main__":
    sequence_file = "data/raw/covid_raw_sequence.fasta"
    metadata_file = "data/raw/covid_meta_data.tsv"
    lineage_file = "data/raw/covid_variant.xlsx"

    covid = CovidFastaProcessor(sequence_file, metadata_file, lineage_file, output_dir="data/clean")
    covid.process("split", "statistic")
