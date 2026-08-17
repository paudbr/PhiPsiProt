#!/usr/bin/env python3

import argparse


def get_chain_residues(pdb_path, chain):
    residues = set()

    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith("ATOM") and line[21] == chain:
                residues.add(int(line[22:26]))

    if not residues:
        raise ValueError(f"No residues found for chain {chain}")

    return sorted(residues)


def continuous_blocks(residues):
    blocks = []
    start = residues[0]
    prev = residues[0]

    for r in residues[1:]:
        if r == prev + 1:
            prev = r
        else:
            blocks.append((start, prev))
            start = r
            prev = r

    blocks.append((start, prev))
    return blocks


def build_region_redesign_contig(pdb_path, chain, redesign_region):
    residues = get_chain_residues(pdb_path, chain)
    blocks = continuous_blocks(residues)

    redesign_start, redesign_end = map(int, redesign_region.split("-"))
    redesign_len = redesign_end - redesign_start + 1

    contig_parts = []

    for start, end in blocks:

        if end < redesign_start or start > redesign_end:
            contig_parts.append(f"{chain}{start}-{end}")
            continue

        if start < redesign_start:
            contig_parts.append(f"{chain}{start}-{redesign_start - 1}")

        contig_parts.append(f"{redesign_len}-{redesign_len}")

        if end > redesign_end:
            contig_parts.append(f"{chain}{redesign_end + 1}-{end}")

    return "/".join(contig_parts)


def build_full_contig(pdb_path, chain):
    residues = get_chain_residues(pdb_path, chain)
    blocks = continuous_blocks(residues)

    return "/".join(
        f"{chain}{start}-{end}"
        for start, end in blocks
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--chain", default="A")
    parser.add_argument(
        "--design_strategy",
        required=True,
        choices=["partial_diffusion", "local_redesign"],
    )
    parser.add_argument("--redesign_region", default=None)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.design_strategy == "local_redesign":
        if not args.redesign_region:
            raise ValueError("--redesign_region is required for local_redesign")

        contig = build_region_redesign_contig(
            args.pdb,
            args.chain,
            args.redesign_region,
        )

    elif args.design_strategy == "partial_diffusion":
        if args.redesign_region:
            contig = build_region_redesign_contig(
                args.pdb,
                args.chain,
                args.redesign_region,
            )
        else:
            contig = build_full_contig(args.pdb, args.chain)

    with open(args.output, "w") as out:
        out.write(contig + "\n")

    print(contig)


if __name__ == "__main__":
    main()