# oncoprint-otp-project
a set of scripts to prepare a project organised by OTP for oncoprint plotting


# How to run

```
echo "modify the perl script project path, preffix, suffix and subfolder paths"
perl make_oncoprint_table.pl -o MY_OUT_PREFFIX
module loads R/3.3.1
R -f make_oncoprint_plot.R --no-save --no-restore --args MY_OUT_PREFFIX.VERSION.min6.kataegis6.sv100000.cnv0.3.onco_print.tsv MY_OUT_PREFFIX.VERSION.min6.kataegis6.sv100000.cnv0.3.sample_info.tsv
```


## Prerequisites

- A project directory structure setup by https://github.com/naveedishaque/otp-project-softlinker
- Developed using perl 5, version 26, subversion 1 (v5.26.1) built for x86_64-linux-gnu-thread-multi
- Conda environment with prequisite tools

## conda setup

[TODO] `conda install perl=5.26.2 perl-datetime=1.42 perl-getopt-long=2.50 bedtools=2.16.2 r=3.4.0 bioconductor-yapsa`

## Usage

```
  perl make_oncoprint_table.pl -i [input directory] -o [output prefix]

```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

 - v1.0.0: first working version (equivalent to v0.7)

## Authors

Naveed Ishaque

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
../otp-project-softlinker/README.md (END)                                                                                                                                                         
