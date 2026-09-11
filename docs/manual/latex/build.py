#!/usr/bin/env python3
"""Build the authoritative OCRT LaTeX manuals. No Markdown conversion is performed."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--language',choices=['ko','en','all'],default='all')
    p.add_argument('--out-dir',type=Path,help='PDF destination; default: docs/manual')
    p.add_argument('--wsl-distro',default='Ubuntu',help='Windows fallback distribution if xelatex is absent')
    args=p.parse_args()
    root=Path(__file__).resolve().parent
    target=(args.out_dir or root.parent).resolve()
    target.mkdir(parents=True,exist_ok=True)
    native=shutil.which('xelatex')
    use_wsl=not native and os.name=='nt' and shutil.which('wsl')
    if not native and not use_wsl:
        p.error('XeLaTeX is required. See README.md for TeX and font dependencies.')
    if use_wsl:
        unix_root=subprocess.check_output(['wsl','-d',args.wsl_distro,'--exec','wslpath','-a','-u',str(root)],text=True).strip()
        prefix=['wsl','-d',args.wsl_distro,'--cd',unix_root,'--exec','xelatex']
    else:
        prefix=[native]
    for language in (['ko','en'] if args.language=='all' else [args.language]):
        name=f'ocrt_manual_{language}'
        build=root/'build'/language
        build.mkdir(parents=True,exist_ok=True)
        command=prefix+['-interaction=nonstopmode','-halt-on-error','-file-line-error',
                        f'-output-directory=build/{language}',name+'.tex']
        for pass_number in range(1,4):
            print(f'{language.upper()}: XeLaTeX pass {pass_number}/3',flush=True)
            result=subprocess.run(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                  text=True,encoding='utf-8',errors='replace')
            (build/f'pass{pass_number}.txt').write_text(result.stdout,encoding='utf-8')
            if result.returncode:
                print(result.stdout[-6000:],file=sys.stderr)
                raise SystemExit(f'Build failed; inspect {build}')
        log=(build/(name+'.log')).read_text(encoding='utf-8',errors='replace')
        error_markers=('Missing character:', 'Undefined control sequence',
                       'ignored error:', 'There were undefined references',
                       'There were multiply-defined labels', 'Rerun to get cross-references right')
        problems=[line for line in log.splitlines() if any(marker in line for marker in error_markers)]
        if problems:
            raise SystemExit('Typesetting or reference errors:\n'+'\n'.join(problems[:25]))
        pdf=build/(name+'.pdf')
        if not pdf.exists(): raise SystemExit('XeLaTeX did not produce the expected PDF')
        destination=target/f'OCRT_Manual_{language.upper()}.pdf'
        shutil.copy2(pdf,destination)
        warnings=sum('Overfull \\hbox' in line or 'Overfull \\vbox' in line for line in log.splitlines())
        print(f'Written: {destination} ({warnings} overfull-box warnings; inspect layout before publication)',flush=True)


if __name__=='__main__':
    main()
