# On a workstation or a cluster

The same pipeline runs on a workstation and on a SLURM cluster: only the
profile and a small config file change. Yeast-sized genomes are small, so most
steps fit anywhere; the few that do not are named below.

## On a workstation

Docker (or Apptainer) and the local executor. Tell Nextflow how much of the
machine it may take, in a config file of your own:

```groovy title="workstation.config"
process {
    resourceLimits = [ cpus: 32, memory: 128.GB, time: 72.h ]
}
executor {
    cpus   = 32
    memory = 128.GB
}
```

```bash
nextflow run js2264/plexembler -profile docker -c workstation.config \
    --input species.csv --libraries libraries.csv --outdir results
```

With an NVIDIA GPU and the NVIDIA container toolkit, add `gpu` to the profile:
the steps labelled `process_gpu` (instaGRAAL's MCMC mode, Helixer,
`dorado polish`) get it (`--gpus all` for Docker, `--nv` for Apptainer).

## On a SLURM cluster

Apptainer (or Singularity), and the SLURM executor. If your institution has a
profile in [nf-core/configs](https://nf-co.re/configs/), use it: it sets the
executor, the queues and the container cache. On Institut Pasteur's Maestro:

```bash
nextflow run js2264/plexembler -profile maestro,normal \
    --input species.csv --libraries libraries.csv --outdir results
```

Otherwise, a config file of your own:

```groovy title="slurm.config"
process {
    executor       = 'slurm'
    queue          = 'common'
    resourceLimits = [ cpus: 64, memory: 512.GB, time: 7.d ]
    withLabel: 'process_gpu' {
        queue          = 'gpu'
        clusterOptions = '--gres=gpu:1'
    }
}
executor {
    queueSize       = 100
    submitRateLimit = '10/1min'
}
apptainer {
    enabled  = true
    cacheDir = '/shared/containers/apptainer'
}
```

```bash
nextflow run js2264/plexembler -profile apptainer,gpu -c slurm.config \
    --input species.csv --libraries libraries.csv --outdir results
```

Run Nextflow itself in `tmux` or `screen`, or as a SLURM job of its own: it
waits on the jobs it submits.

## What needs more than the rest

| Step | What it needs | When it runs |
|---|---|---|
| FCS-GX (foreign DNA) | its database, about 500 GB, best held in memory | when `--fcsgx_db` is given |
| InterProScan | its database, about 100 GB | when `--interproscan_db` is given |
| instaGRAAL, `--instagraal_mode default` | an NVIDIA GPU | `cpu` mode needs none |
| Helixer, `dorado polish` | a GPU, for speed | optional steps |

These come with their milestones ([implementation plan](../plan.md)). Databases
are given as paths, downloaded once to shared storage; the small ones
(BUSCO's lineages, OrthoDB's proteins, Rfam) are downloaded by the pipeline
into a cache folder when no path is given.

## Resuming and cleaning up

`-resume` reuses every finished step. The work directory (`-w`, by default
`./work`) holds every step's files: put it on scratch storage on a cluster,
and delete it once the results are checked.
