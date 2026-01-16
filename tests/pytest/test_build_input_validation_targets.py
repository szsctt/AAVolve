from aavolve.snakemake_helpers import build_group_report_targets, build_input_validation_targets


def test_build_input_validation_targets_deduplicates(samples_df):
    input_validation_map, input_validation_samples, input_validation_targets = build_input_validation_targets(
        samples_df
    )

    # At least one target, and it matches map size.
    assert len(input_validation_targets) == len(input_validation_map)

    # Ensure each target corresponds to a pair_id.
    for target in input_validation_targets:
        assert target.startswith("out/qc/input-checks/")
        assert target.endswith(".txt")
        pair_id = target.split("/")[-1].replace(".txt", "")
        assert pair_id in input_validation_map
        assert pair_id in input_validation_samples


def test_build_group_report_targets(samples_df):
    input_validation_map, input_validation_samples, group_report_targets = build_group_report_targets(samples_df)

    assert len(group_report_targets) == len(input_validation_map)
    for target in group_report_targets:
        assert target.startswith("out/qc/group_reports/")
        assert target.endswith("_report.html")
        pair_id = target.split("/")[-1].replace("_report.html", "")
        assert pair_id in input_validation_map
        assert pair_id in input_validation_samples
