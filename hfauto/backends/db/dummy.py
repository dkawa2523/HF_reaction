from __future__ import annotations


class DummyDBProvider:
    name = "dummy"

    def __init__(self, **kwargs):
        self.config = kwargs

    def enrich(self, molecule_artifact):
        data = molecule_artifact.data.copy()
        smiles = data.get("canonical_smiles", "")
        data.setdefault("identity", {})
        data["identity"].update(
            {
                "provider": "dummy",
                "identity_confidence": "unknown" if not smiles else "structure_only",
                "identity_conflict": False,
            }
        )
        data.setdefault("public_data", {})
        data["public_data"].setdefault("gas_process", {})
        data["public_data"]["gas_process"].update(
            {
                "gas_process_feasibility": "not_evaluated_dummy",
                "ehs_review_flag": "not_evaluated_dummy",
            }
        )
        data.setdefault("db_provider_status", {})[self.name] = {"status": "success"}
        return data
