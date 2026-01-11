import yaml
import subprocess
import datetime
import platform

class RunManifest:
    def __init__(self, output_path="results/run.yaml"):
        self.output_path = output_path
        self.data = {
            "timestamp": datetime.datetime.now().isoformat(),
            "environment": {
                "os": platform.system(),
                "python_version": platform.python_version()
            },
            "tools": {},
            "parameters": {}
        }

    def capture_tool_versions(self):
        """Records versions of external bioinformatics tools[cite: 20, 177]."""
        # Capture MAFFT version
        try:
            mafft_v = subprocess.check_output(["mafft", "--version"], stderr=subprocess.STDOUT).decode()
            self.data["tools"]["mafft"] = mafft_v.strip()
        except:
            self.data["tools"]["mafft"] = "not found"

        # Capture DSSP version
        try:
            dssp_v = subprocess.check_output(["dssp", "--version"], stderr=subprocess.STDOUT).decode()
            self.data["tools"]["dssp"] = dssp_v.strip()
        except:
            self.data["tools"]["dssp"] = "not found"

    def add_parameter(self, stage, params):
        """Logs parameters used in specific pipeline stages[cite: 178]."""
        self.data["parameters"][stage] = params

    def save(self):
        with open(self.output_path, 'w') as f:
            yaml.dump(self.data, f, default_flow_style=False)
        print(f"📄 Run manifest saved to {self.output_path}")