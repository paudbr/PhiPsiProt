from pathlib import Path
import re
import os

file_path = Path(
    "/chai1_venv/lib/python3.10/site-packages/chai_lab/data/dataset/embeddings/esm.py"
)

text = file_path.read_text()

# asegurar import os
if "import os" not in text:
    text = "import os\n" + text

pattern = r"with esm_model\(\s*device\s*=\s*device\s*\)\s*as\s*model:"

replacement = """with esm_model(
        device=(
            device
            if not (device.type == "cuda" and os.environ.get("CHAI_FORCE_ESM_CPU", "0") == "1")
            else torch.device("cpu")
        )
    ) as model:"""

if not re.search(pattern, text):
    raise RuntimeError("PATCH FAILED: esm_model(device=device) not found")

text = re.sub(pattern, replacement, text)

file_path.write_text(text)

print("OK: esm_model vGPU-safe patch applied")