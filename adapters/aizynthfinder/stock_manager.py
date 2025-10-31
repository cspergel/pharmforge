#!/usr/bin/env python3
"""
Stock Database Manager for AiZynthFinder

Allows users to easily create and manage custom stock databases
for retrosynthesis planning.

Features:
- Create stock databases from simple CSV files
- Add/remove compounds interactively
- Validate SMILES strings
- Generate HDF5 stock files compatible with AiZynthFinder
"""

import os
import logging
from typing import List, Dict, Set, Optional
from pathlib import Path

logger = logging.getLogger(__name__)


class StockManager:
    """
    Manages custom stock databases for AiZynthFinder.

    Allows users to define available compounds via CSV and generates
    HDF5 stock files that AiZynthFinder can use.
    """

    # Common lab reagents (expandable by users)
    DEFAULT_LAB_STOCK = {
        # Solvents
        "O": "water",
        "CO": "methanol",
        "CCO": "ethanol",
        "CC(C)O": "isopropanol",
        "CC(=O)C": "acetone",
        "CCOC(C)=O": "ethyl acetate",
        "ClCCl": "dichloromethane",
        "ClC(Cl)Cl": "chloroform",
        "c1ccccc1": "benzene",
        "Cc1ccccc1": "toluene",
        "CC#N": "acetonitrile",
        "CN(C)C=O": "DMF",
        "CS(C)=O": "DMSO",
        "C1CCOC1": "THF",
        "C1CCCCC1": "cyclohexane",
        "CCCCCC": "hexane",

        # Acids
        "CC(=O)O": "acetic acid",
        "O=C(O)C(O)C(O)C(O)C(O)CO": "gluconic acid",
        "ClC(Cl)(Cl)C(=O)O": "trichloroacetic acid",
        "O=C(O)C(F)(F)F": "trifluoroacetic acid",
        "O=S(=O)(O)O": "sulfuric acid",
        "Cl": "hydrochloric acid",
        "O=C(O)c1ccccc1": "benzoic acid",

        # Bases
        "C(C)(C)(C)O[K]": "potassium tert-butoxide",
        "[Na+].[OH-]": "sodium hydroxide",
        "[K+].[OH-]": "potassium hydroxide",
        "CCN(CC)CC": "triethylamine",
        "C1CCNCC1": "piperidine",
        "c1cccnc1": "pyridine",

        # Common reagents
        "CC(=O)OC(=O)C": "acetic anhydride",
        "O=C(Cl)C(=O)Cl": "oxalyl chloride",
        "ClC(Cl)=O": "phosgene",
        "O=C=O": "carbon dioxide",
        "ClCCCl": "1,2-dichloroethane",

        # Building blocks - benzene derivatives
        "Nc1ccccc1": "aniline",
        "Oc1ccccc1": "phenol",
        "c1ccc(C=O)cc1": "benzaldehyde",
        "c1ccc(C(=O)Cl)cc1": "benzoyl chloride",
        "c1ccc(CCl)cc1": "benzyl chloride",
        "c1ccc(CO)cc1": "benzyl alcohol",
        "c1ccc(C(=O)O)cc1": "benzoic acid",
        "c1ccc(N)cc1": "aniline",

        # Building blocks - aliphatic
        "C=C": "ethylene",
        "C=CC": "propylene",
        "C=O": "formaldehyde",
        "CC=O": "acetaldehyde",
        "CC(=O)C": "acetone",
        "CCC=O": "propionaldehyde",
        "CCCC": "butane",
        "C=CCC": "1-butene",

        # Halides
        "CCCBr": "1-bromopropane",
        "CCBr": "bromoethane",
        "CBr": "bromomethane",
        "CCI": "iodomethane",
        "CCCl": "chloroethane",

        # Amines
        "CN": "methylamine",
        "CCN": "ethylamine",
        "CCCN": "propylamine",
        "CNC": "dimethylamine",
        "CCNCC": "diethylamine",

        # Alcohols
        "C": "methane",
        "CC": "ethane",
        "CCC": "propane",
        "CCCC": "butane",
        "CCCCC": "pentane",
        "CCCCCC": "hexane",
    }

    def __init__(self, data_dir: str = "/app/aizynthfinder_data"):
        """
        Initialize stock manager.

        Args:
            data_dir: Directory for stock database files
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # Paths
        self.stock_csv = self.data_dir / "custom_stock.csv"
        self.stock_hdf5 = self.data_dir / "custom_stock.hdf5"

    def validate_smiles(self, smiles: str) -> bool:
        """
        Validate SMILES string.

        Args:
            smiles: SMILES string

        Returns:
            True if valid, False otherwise
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation failed for '{smiles}': {e}")
            return False

    def canonicalize_smiles(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to canonical form.

        Args:
            smiles: SMILES string

        Returns:
            Canonical SMILES or None if invalid
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
            return None
        except Exception:
            return None

    def load_stock_from_csv(self, csv_path: Optional[str] = None) -> Dict[str, str]:
        """
        Load stock compounds from CSV file.

        CSV format: smiles,name (one compound per line)

        Args:
            csv_path: Path to CSV file (default: custom_stock.csv)

        Returns:
            Dictionary mapping SMILES to compound names
        """
        if csv_path is None:
            csv_path = self.stock_csv

        stock = {}

        if not os.path.exists(csv_path):
            logger.info(f"No CSV file found at {csv_path}, using defaults")
            return self.DEFAULT_LAB_STOCK.copy()

        try:
            with open(csv_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue

                    # Parse CSV
                    parts = line.split(',')
                    if len(parts) < 1:
                        logger.warning(f"Line {line_num}: Invalid format")
                        continue

                    smiles = parts[0].strip()
                    name = parts[1].strip() if len(parts) > 1 else "unnamed"

                    # Validate SMILES
                    canonical = self.canonicalize_smiles(smiles)
                    if canonical:
                        stock[canonical] = name
                    else:
                        logger.warning(f"Line {line_num}: Invalid SMILES '{smiles}'")

            logger.info(f"Loaded {len(stock)} compounds from {csv_path}")

        except Exception as e:
            logger.error(f"Failed to load CSV: {e}")
            return self.DEFAULT_LAB_STOCK.copy()

        return stock

    def save_stock_to_csv(self, stock: Dict[str, str], csv_path: Optional[str] = None):
        """
        Save stock compounds to CSV file.

        Args:
            stock: Dictionary mapping SMILES to names
            csv_path: Path to CSV file (default: custom_stock.csv)
        """
        if csv_path is None:
            csv_path = self.stock_csv

        try:
            with open(csv_path, 'w') as f:
                f.write("# Custom Stock Database for AiZynthFinder\n")
                f.write("# Format: smiles,name\n")
                f.write("# Edit this file to add/remove available compounds\n\n")

                for smiles, name in sorted(stock.items(), key=lambda x: x[1]):
                    f.write(f"{smiles},{name}\n")

            logger.info(f"Saved {len(stock)} compounds to {csv_path}")

        except Exception as e:
            logger.error(f"Failed to save CSV: {e}")
            raise

    def create_hdf5_stock(self, stock: Dict[str, str], hdf5_path: Optional[str] = None) -> bool:
        """
        Create HDF5 stock database from SMILES dictionary.

        Args:
            stock: Dictionary mapping SMILES to names
            hdf5_path: Path to output HDF5 file (default: custom_stock.hdf5)

        Returns:
            True if successful, False otherwise
        """
        if hdf5_path is None:
            hdf5_path = self.stock_hdf5

        try:
            from rdkit import Chem
            import pandas as pd

            # Convert to list of InchiKey entries
            inchi_keys = []
            for smiles in stock.keys():
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    inchi_key = Chem.MolToInchiKey(mol)
                    inchi_keys.append(inchi_key)

            # Create DataFrame with inchi_key column (same format as ZINC stock)
            df = pd.DataFrame({'inchi_key': inchi_keys})

            # Save as HDF5 with 'table' key (AiZynthFinder format)
            df.to_hdf(hdf5_path, key='table', mode='w', format='table')

            logger.info(f"Created HDF5 stock with {len(inchi_keys)} compounds at {hdf5_path}")
            return True

        except ImportError as e:
            logger.error(f"Required library not available: {e}")
            return False
        except Exception as e:
            logger.error(f"Failed to create HDF5 stock: {e}")
            return False

    def initialize_default_stock(self) -> bool:
        """
        Initialize default lab stock database.

        Creates both CSV and HDF5 files with common lab reagents.

        Returns:
            True if successful, False otherwise
        """
        logger.info("Initializing default lab stock database...")

        # Save default stock to CSV
        self.save_stock_to_csv(self.DEFAULT_LAB_STOCK)

        # Create HDF5 from defaults
        success = self.create_hdf5_stock(self.DEFAULT_LAB_STOCK)

        if success:
            logger.info(f"✓ Default stock initialized with {len(self.DEFAULT_LAB_STOCK)} compounds")
            logger.info(f"  CSV: {self.stock_csv}")
            logger.info(f"  HDF5: {self.stock_hdf5}")

        return success

    def update_stock_database(self, csv_path: Optional[str] = None) -> bool:
        """
        Update HDF5 stock database from CSV file.

        Args:
            csv_path: Path to CSV file (default: custom_stock.csv)

        Returns:
            True if successful, False otherwise
        """
        logger.info("Updating stock database from CSV...")

        # Load from CSV
        stock = self.load_stock_from_csv(csv_path)

        # Create HDF5
        success = self.create_hdf5_stock(stock)

        if success:
            logger.info(f"✓ Stock database updated with {len(stock)} compounds")

        return success

    def add_compound(self, smiles: str, name: str) -> bool:
        """
        Add a compound to the stock database.

        Args:
            smiles: SMILES string
            name: Compound name

        Returns:
            True if successful, False otherwise
        """
        # Validate and canonicalize
        canonical = self.canonicalize_smiles(smiles)
        if not canonical:
            logger.error(f"Invalid SMILES: {smiles}")
            return False

        # Load current stock
        stock = self.load_stock_from_csv()

        # Add compound
        stock[canonical] = name

        # Save updated stock
        self.save_stock_to_csv(stock)

        # Regenerate HDF5
        return self.create_hdf5_stock(stock)

    def remove_compound(self, smiles: str) -> bool:
        """
        Remove a compound from the stock database.

        Args:
            smiles: SMILES string

        Returns:
            True if successful, False otherwise
        """
        # Canonicalize
        canonical = self.canonicalize_smiles(smiles)
        if not canonical:
            logger.error(f"Invalid SMILES: {smiles}")
            return False

        # Load current stock
        stock = self.load_stock_from_csv()

        # Remove compound
        if canonical in stock:
            del stock[canonical]
            logger.info(f"Removed {canonical}")
        else:
            logger.warning(f"Compound not found: {canonical}")
            return False

        # Save updated stock
        self.save_stock_to_csv(stock)

        # Regenerate HDF5
        return self.create_hdf5_stock(stock)

    def get_stock_info(self) -> Dict[str, any]:
        """
        Get information about the current stock database.

        Returns:
            Dictionary with stock statistics
        """
        stock = self.load_stock_from_csv()

        return {
            "total_compounds": len(stock),
            "csv_path": str(self.stock_csv),
            "hdf5_path": str(self.stock_hdf5),
            "csv_exists": self.stock_csv.exists(),
            "hdf5_exists": self.stock_hdf5.exists(),
            "compounds": stock
        }


def main():
    """CLI interface for stock manager"""
    import argparse

    parser = argparse.ArgumentParser(description="Manage AiZynthFinder stock databases")
    parser.add_argument('--init', action='store_true', help='Initialize default stock database')
    parser.add_argument('--update', action='store_true', help='Update HDF5 from CSV')
    parser.add_argument('--add', nargs=2, metavar=('SMILES', 'NAME'), help='Add compound')
    parser.add_argument('--remove', metavar='SMILES', help='Remove compound')
    parser.add_argument('--info', action='store_true', help='Show stock info')
    parser.add_argument('--data-dir', default='/app/aizynthfinder_data', help='Data directory')

    args = parser.parse_args()

    manager = StockManager(data_dir=args.data_dir)

    if args.init:
        manager.initialize_default_stock()
    elif args.update:
        manager.update_stock_database()
    elif args.add:
        smiles, name = args.add
        manager.add_compound(smiles, name)
    elif args.remove:
        manager.remove_compound(args.remove)
    elif args.info:
        info = manager.get_stock_info()
        print(f"\nStock Database Info:")
        print(f"  Total compounds: {info['total_compounds']}")
        print(f"  CSV file: {info['csv_path']} ({'exists' if info['csv_exists'] else 'missing'})")
        print(f"  HDF5 file: {info['hdf5_path']} ({'exists' if info['hdf5_exists'] else 'missing'})")
        print(f"\nAvailable compounds:")
        for smiles, name in sorted(info['compounds'].items(), key=lambda x: x[1])[:20]:
            print(f"    {name:<30} {smiles}")
        if info['total_compounds'] > 20:
            print(f"    ... and {info['total_compounds'] - 20} more")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
