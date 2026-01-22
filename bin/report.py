#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
import pysam

from typing import List, Type
import base64

from bokeh.models import Title
from dominate.tags import (
    h1, p, div, section, html_tag,
    a, button, img, ul, li, strong,
    span
)

import ezcharts as ezc

from ezcharts.layout.snippets import Tabs, Grid
from ezcharts.layout.snippets.table import DataTable
from ezcharts.layout.snippets.section import Section
from ezcharts.layout.snippets.tabs import ITabsClasses, ITabsStyles
from ezcharts.layout.snippets.banner import (
    IBannerStyles, IBannerClasses
)

from ezcharts.layout.base import Snippet
from ezcharts.layout.util import cls, css

from ezcharts.components.reports import Report, labs
from ezcharts.components.ezchart import EZChart
from ezcharts.components.theme import LAB_body_resources, LAB_head_resources


BACKGROUND_COLOR = '#1d1d1d'
PLOT_COLOR = None


class NoMarginTabsClasses(ITabsClasses):
    """Override tab classes to remove margins."""
    tab_buttons_list: str = cls("nav", "nav-tabs")  # Removed "mb-2"
    tab_contents_container: str = cls("tab-content", "p-0")  # Set padding to 0


class NoMarginTabsStyles(ITabsStyles):
    """Override tab styles to remove margins."""
    tab_button: str = css(
        "margin-bottom: 0",  # Remove negative margin
        "font-weight: 600",
        "cursor: pointer",
        "border-color: transparent"
    )
    tab_button_active: str = css(
        "border-bottom: 2px solid #0079a4",
        "color: #0079a4!important"
    )


class NoMarginTabs(Tabs):
    """Custom Tabs class without extra margins."""
    def __init__(self) -> None:
        super().__init__(styles=NoMarginTabsStyles(), classes=NoMarginTabsClasses())


class OprBanner(Snippet):
    """A styled div tag containing a heading and badges."""

    TAG = 'div'

    def __init__(self, report_title: str, pipeline_version: str) -> None:
        # TODO: replace the banner styles and classes?
        styles: IBannerStyles = IBannerStyles()
        classes: IBannerClasses = IBannerClasses()
        # Override the default background color
        styles.container = css(
            f"background-color: {BACKGROUND_COLOR} !important;",
            "padding: 15px;",
            "border-radius: 5px;"
        )

        super().__init__(
            styles=styles,
            classes=classes,
            className=classes.container,
            style=styles.container
        )
        with self:
            with div(className=self.classes.inner, style=self.styles.inner):
                h1(report_title)
                p(f'Pipeline version {pipeline_version}')
                p('Research use only')


class OprLogo(div):
    """OPR logo element."""
    def __init__(self) -> None:
        fname_logo = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'OPR_logo_v01-light.cropped.png'
        )

        # Convert image to Base64
        with open(fname_logo, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode("utf-8")

        super().__init__(
            img(src=f"data:image/png;base64,{encoded_string}", style="height: 75px;", alt="OPR Logo"),
            tagname='div',
            className="d-flex"
        )


class OprNavigation(Snippet):
    """A styled nav component for use in a Report."""

    TAG = 'nav'

    def __init__(
        self,
        logo: Type[html_tag],
        groups: List[str],
        header_height: int = 75,
        classes: labs.ILabsNavigationClasses = labs.ILabsNavigationClasses()
    ) -> None:
        spacer = div(
            className=classes.spacer,
            style=f"margin-top: {header_height}px;"
        )
        super().__init__(
            styles=None,
            classes=classes,
            style=f"min-height: {header_height}px; background-color: {BACKGROUND_COLOR} !important;",
            className=classes.container
        )
        spacer.add(self)
        with self:
            with div(className=self.classes.inner):
                with a(href="https://github.com/OPR-group-BNITM",
                       className=self.classes.logo):
                    logo()
                button(
                    "Jump to section... ",
                    cls=self.classes.dropdown_btn,
                    type="button",
                    id="dropdownMenuButton",
                    data_bs_toggle="dropdown",
                    aria_haspopup="true",
                    aria_expanded="false"
                )
                ngroups = len(groups)
                with div(className=self.classes.dropdown_menu):
                    for count, group in enumerate(groups):
                        setattr(
                            self, group,
                            div(className='', __pretty=False)
                        )
                        if count != ngroups - 1:
                            div(cls="dropdown-divider")

    def add_link(
        self,
        group: str,
        link_title: str,
        link_href: str
    ) -> None:
        """Add a header nav link to the header links list."""
        group_list = getattr(self, group)
        with group_list:
            a(
                link_title,
                href=link_href,
                className=self.classes.dropdown_item_link
            )


class OprReport(Report):
    """A basic OPR-themed report."""

    def __init__(self, title, pipeline_version):
        super().__init__(
            report_title=title,
            head_resources=LAB_head_resources,
            body_resources=LAB_body_resources
        )
        with self.header:
            self.header.attributes["style"] = f"background-color: {BACKGROUND_COLOR} !important;"
            self.nav = OprNavigation(logo=OprLogo, groups=['main', 'meta'])
            self.intro_content = section(
                id="intro-content",
                role="region",
                style=f"background-color: {BACKGROUND_COLOR}  !important;"
            )
            with self.intro_content:
                self.banner = OprBanner(title, pipeline_version)

        with self.main:
            self.main_content = section(id="main-content", role="region")

    def add_section(
        self,
        title: str,
        link: str,
        overflow: bool = False
    ) -> Section:
        """Add a section to the main_content region."""
        href = self.get_uid('Section')
        self.nav.add_link('main', link, f'#{href}')
        with self.main_content:
            return Section(href, title, overflow=overflow)


def histogram_plot(data, title, xaxis_label):

    nbins = 300
    min_val = data.min()
    max_val = data.max()
    bin_width = (max_val - min_val) / nbins
    mean = data.mean()
    median = data.median()

    plot = ezc.histplot(
        data=data,
        binwidth=bin_width,
        binrange=(min_val, max_val),
        color=PLOT_COLOR
    )
    subtitle = f"Mean: {mean:.1f}. Median: {median:.1f}"

    plot._fig.xaxis.axis_label = xaxis_label
    plot._fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), 'above')
    plot._fig.add_layout(Title(text=title, text_font_size="1.5em"), 'above')
    
    plot._fig.x_range.start = min_val
    plot._fig.x_range.end = max_val

    plot._fig.yaxis.axis_label = 'Number of reads'

    return plot


class Legends:

    def __init__(self, descriptions: dict, default_title: str = None):
        self.descriptions = descriptions
        self.default_title = default_title if default_title else ""

    def legend(self, columns: list, title: str = None):
        with div(style="font-size: 0.9em; margin-top: 1em; color: #5c5c5c;"):
            title = title if title else self.default_title
            p(title, style="font-size: 1em; font-weight: bold; margin-bottom: 0.3em;")
            with ul(style="padding-left: 1.2em; margin: 0;"):
                for col in columns:
                    with li():
                        strong(f"{col}: ")
                        span(self.descriptions[col])


def genome_properties(fname):
    genome = SeqIO.read(fname, 'fasta')
    len_genome = len(genome.seq)
    recoverage = 1 - str(genome.seq).count('N') / len_genome
    properties = {
        'Genome Length': len_genome,
        'Recovered': recoverage,
    }
    return pd.DataFrame([properties])


def vcf_header_columns(vcf, header_name):
    return [
        head.strip().strip("'")
        for head in vcf.header.info[header_name].description.split(':')[1].split('|')
    ]


def getval(val):
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def vcf_info(variant, header_name, headers):
    if header_name in variant.info:
        vals = [
            getval(val.strip())
            for val in variant.info[header_name][0].split('|')
        ]
        return dict(zip(headers, vals))
    return {}


def read_variants(fname):
    lines = []
    with pysam.VariantFile(fname, 'r') as vcf:
        cols_ann = vcf_header_columns(vcf, 'ANN')
        cols_basecounts = vcf_header_columns(vcf, 'BASECOUNTS')
        for variant in vcf:
            variant_info = {
                'Reference': variant.chrom,
                'Position': variant.pos,
                'Ref': variant.ref,
                'Alt': ','.join(str(alt) for alt in variant.alts),
                'Qual': variant.qual,
                'Filter': variant.filter.keys()
            }
            annotations = vcf_info(variant, 'ANN', cols_ann)
            basecounts = vcf_info(variant, 'BASECOUNTS', cols_basecounts)
            if variant_info['Alt'].strip() != 'N':
                lines.append(variant_info | annotations | basecounts)
    return pd.DataFrame(lines)


def html_report(
        pipeline_version,
        samplename,
        genome_props,
        variants,
        fname_out 
):
    report = OprReport(
        title=f"MOPX report for sample {samplename}",
        pipeline_version=pipeline_version
    )

    with report.add_section("Genome properties", "Genome properties"):
        cols = ['Genome Length', 'Recovered']
        DataTable.from_pandas(
            genome_props[cols],
            use_index=False,
            export=False
        )
        Legends(
            {
                'Genome Length': 'Length of the output genome',
                'Recovered': 'Amount of the genome recovered in % that are not N.',
            },
            'Genome recovery'
        ).legend(cols)

    with report.add_section("Variants", "Variants"):
        cols = [
            'Position',
            'Ref',
            'Alt',
            'Annotation',
            'Gene_Name',
            'HGVS.p',
        ]
        DataTable.from_pandas(
            variants[cols],
            use_index=False,
            export=False
        )
        Legends(
            {
                'Position': 'Genomic position of the variant.',
                'Ref': 'Reference that is mutated.',
                'Alt': 'Mutated seqeunce.',
                'Annotation': 'Class of the mutation impact.',
                'Gene_Name': 'Name of the gene.',
                'HGVS.p': 'Amino acid sequence change',
            },
            'Variants'
        ).legend(cols)

    report.write(fname_out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samplename', required=True, help='Name of the sample to report.')
    parser.add_argument('--mopx-version', required=True, help='MOPX version used.')
    parser.add_argument('--draft', required=True, help='Consensus genome created by the pipeline.')
    parser.add_argument('--vcf', required=True, help='VCF file with anotated variants.')
    parser.add_argument('--out', required=True, help='Name of the output html report file.')
    args = parser.parse_args()

    genome_props = genome_properties(args.draft)
    variants = read_variants(args.vcf)

    html_report(
        args.mopx_version,
        args.samplename,
        genome_props,
        variants,
        args.out
    )


if __name__ == '__main__':
    main()


# TODO:
# - Report if no bases called!
# - Report if no variants
