% ==========================================================
% METADATA
% ==========================================================
%%%% Template for ULB report %%%%
% By Boris Coquelet %
%
%Dans le cadre du poste aide à la réussite au Bureau des étudiants de Polytechnique%
%2019

% ==========================================================
% CLASS DECLARATION
% ==========================================================
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ULBreportEnhanced}[2025/07/10 v2.3 - Enhanced with Full Spacing Control]
\LoadClass[a4paper,11pt]{report}

% ==========================================================
% PACKAGE LOADING
% ==========================================================

% --- Core & Language ---
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{csquotes}
\RequirePackage[english]{babel}

% --- Layout & Formatting ---
\RequirePackage[top=2.5cm,bottom=2.5cm,left=2.5cm,right=2.5cm]{geometry}
\RequirePackage{setspace}
\RequirePackage{fancyhdr}
\RequirePackage{titlesec}
%\RequirePackage{newtxtext, newtxmath}
%\RequirePackage{fourier}
\usepackage{lmodern}
\RequirePackage[toc,page]{appendix}
\RequirePackage{indentfirst}

% --- Graphics & Color ---
\RequirePackage{graphicx}
\RequirePackage{xcolor}
\RequirePackage{colortbl}
\RequirePackage{tikz}
\RequirePackage{pgfplots}
\RequirePackage{chemfig}
\RequirePackage{subcaption}
\RequirePackage{float}

% --- Math ---
\RequirePackage{mathtools}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{mathrsfs}
\RequirePackage[makeroom]{cancel}
\RequirePackage{esint}

% --- Bibliography ---
\RequirePackage[backend=biber,citestyle=alphabetic,bibstyle=authoryear]{biblatex}

% --- Special Content & Utilities ---
\RequirePackage{tabulaRX}
\RequirePackage{listings}
\RequirePackage{etoolbox}
\RequirePackage[pages=some]{background}
\RequirePackage[most]{tcolorbox}
\RequirePackage{hyperref}
\RequirePackage{enumitem}
\RequirePackage{xkeyval}
\RequirePackage{ifthen}

% ==========================================================
% TIKZ LIBRARIES
% ==========================================================
\usetikzlibrary{mindmap, arrows, matrix, positioning, 3d, calc}
\pgfplotsset{compat=newest}

% ==========================================================
% GLOBAL CONFIGURATIONS
% ==========================================================

% --- Color Palette ---
\definecolor{primary}{HTML}{005A9C}
\definecolor{secondary}{HTML}{4E5D6C}
\definecolor{accent}{HTML}{E63946}
\definecolor{bgcolor}{HTML}{F1F1F1}

\hypersetup{
    colorlinks=true,
    linkcolor=black,
    filecolor=primary,
    urlcolor=primary,
    pdfpagemode=FullScreen,
}

% --- Typesetting ---
\linespread{1.3}

% --- Page Style (Headers/Footers) ---
\setlength{\headheight}{14pt} % Fix for fancyhdr warning
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textit{\nouppercase{\leftmark}}}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0pt}
\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyfoot[C]{\thepage}%
  \renewcommand{\headrulewidth}{0pt}%
}

% --- Chapter Title Formatting (Full Spacing Control) ---
\titleformat{\chapter}[display]
  {\normalfont\huge\bfseries\color{primary}}
  {\filleft\parbox[b]{2cm}{\fontsize{80}{80}\selectfont\bfseries\color{primary!20}\thechapter}}
  {0pt}
  {\filright}
  [\vspace{1ex}\titlerule]

% Spacing: \titlespacing*{<command>}{<left>}{<before-sep>}{<after-sep>}
% <before-sep> is negative to pull the title up.
% <after-sep> controls space between title and content.
\titlespacing*{\chapter}{0pt}{-10pt}{20pt}


% --- Appendix Naming ---
\renewcommand{\appendixtocname}{Annexes}
\renewcommand{\appendixpagename}{Annexes}

% --- Bibliography Formatting ---
\DeclareFieldFormat{labelalphawidth}{\mkbibbrackets{#1}}
\defbibenvironment{bibliography}
   {\list
      {\printtext[labelalphawidth]{%
         \printfield{prefixnumber}%
         \printfield{labelalpha}%
         \printfield{extraalpha}}}
      {\setlength{\labelwidth}{\labelalphawidth}%
       \setlength{\leftmargin}{\labelwidth}%
       \setlength{\labelsep}{\biblabelsep}%
       \addtolength{\leftmargin}{\labelsep}%
       \setlength{\itemsep}{\bibitemsep}%
       \setlength{\parsep}{\bibparsep}}%
       \renewcommand*{\makelabel}[1]{##1\hss}}
   {\endlist}
   {\item}

% ==========================================================
% CUSTOM ENVIRONMENTS
% ==========================================================

% --- Highlight Box Environment ---
\newtcolorbox{highlightbox}[2][]{
  colback=bgcolor,
  colframe=primary,
  fonttitle=\bfseries,
  coltitle=white,
  colbacktitle=primary,
  title=#2,
  enhanced,
  attach boxed title to top center={yshift=-2mm},
  sharp corners,
  boxrule=1pt,
  {#1}
}

% ==========================================================
% CUSTOM COMMANDS
% ==========================================================
\makeatletter

% --- Background Seal Command (from Original version) ---
\newcommand{\sceau}[1]{\backgroundsetup{
scale=1,
color=black,
opacity=0.03,
angle=10,
position={12cm,-22cm},
contents={%
  \includegraphics[height=20cm,width=20cm,keepaspectratio]{#1}
  }%
}}

% --- Title Page Command (from Original version) ---

% Keys for title page
\define@cmdkey [PRE] {fam} {title}{}
\define@cmdkey [PRE] {fam} {studies}{}
\define@cmdkey [PRE] {fam} {course}{}
\define@cmdkey [PRE] {fam} {author}{}
\define@cmdkey [PRE] {fam} {teacher}{}
\define@cmdkey [PRE] {fam} {date}{}
\define@cmdkey [PRE] {fam} {logo}{}
\define@cmdkey [PRE] {fam} {manyAuthor}[1]{}
% Set default value for the keys (used if not specified)
\presetkeys           [PRE] {fam} {title=titre,
                                   studies=,
                                   course=,
                                   author=Auteurs,
                                   teacher=,
                                   date=,
                                   logo=7,
                                   manyAuthor=8,}{}

% Title command definition
\newcommand*{\titleULB}[1]{\setkeys[PRE]{fam}{#1}
    \begingroup
    \BgThispage
    \thispagestyle{empty}
    % Replace logo with nothing if not specified
    \ifthenelse{\equal{\cmdPRE@fam@logo}{7}}{}{\includegraphics[width=0.87\linewidth,keepaspectratio]{\cmdPRE@fam@logo}\\}
        {\selectfont \cmdPRE@fam@studies}
        
    \begin{center}
    \vspace{1cm}
    \vspace*{\baselineskip}
    {\selectfont\bfseries\cmdPRE@fam@course}\\[0\baselineskip]
    \vspace{2cm}
    \rule{\textwidth}{1.6pt}\vspace*{-\baselineskip}\vspace*{2pt}
    \rule{\textwidth}{0.4pt}\\[\baselineskip]
    \parbox{\textwidth}{\centering\bfseries\selectfont\fontsize{20pt}{24pt}\selectfont \cmdPRE@fam@title}\\[\baselineskip]
    \rule{\textwidth}{0.4pt}\vspace*{-\baselineskip}\vspace{3.2pt}
    \rule{\textwidth}{1.6pt}\\[\baselineskip]
    \vspace{5cm}
    \ifthenelse{\equal{\cmdPRE@fam@manyAuthor}{8}}{
    {\Large\selectfont \cmdPRE@fam@author}\\[\baselineskip]
    {\bfseries\selectfont \cmdPRE@fam@teacher}}{
        \begin{minipage}{0.4\textwidth}
            \begin{flushleft} \large\selectfont
            \cmdPRE@fam@author
            \end{flushleft}
        \end{minipage}
        \begin{minipage}{0.4\textwidth}
            \begin{flushright} \large\selectfont
            \cmdPRE@fam@teacher
            \end{flushright}
        \end{minipage}
    }
    \vfill
    {\selectfont \cmdPRE@fam@date}
    \end{center}
    \vspace*{3\baselineskip}
    \pagenumbering{gobble}
    \newpage
    \tableofcontents
    \newpage
    \pagenumbering{arabic}
\endgroup}

\makeatother
