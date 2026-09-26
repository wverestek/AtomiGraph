# shared test constants
from pathlib import Path

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
PE_DIR = EXAMPLES / "01_PE_chain_basic"
PEEK_DUMP = EXAMPLES / "02_PEEK_one_reaction" / "bonds.reaxff.dump"
PEEK_TYPE_MAP = "1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O"
