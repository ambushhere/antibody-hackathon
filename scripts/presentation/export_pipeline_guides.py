"""Render the Russian and English Markdown guides as searchable, illustrated PDFs."""
import argparse
import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Flowable, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, PageBreak

ROOT = Path(__file__).resolve().parents[2]
INK = colors.HexColor('#10243e')
BLUE = colors.HexColor('#2867b2')
GREEN = colors.HexColor('#267557')
AMBER = colors.HexColor('#a45b08')
W = A4[0] - 96
BASE_URL = 'https://github.com/ambushhere/antibody-hackathon/blob/main/'


def clean(text):
    return text.replace('\u2011', '-').replace('—', '-').replace('–', '-')


def inline(text, source):
    text = html.escape(clean(text))
    def link(match):
        title, target = match.groups()
        if not target.startswith(('http://', 'https://')):
            target = BASE_URL + str((source.parent / target).resolve().relative_to(ROOT))
        return f'<a href="{target}" color="#2867b2">{title}</a>'
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', link, text)
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    return re.sub(r'`([^`]+)`', r'<font color="#285577">\1</font>', text)


class Diagram(Flowable):
    def __init__(self, number, lang, style):
        Flowable.__init__(self)
        self.width, self.height = W, 200
        self.number, self.lang, self.style = number, lang, style

    def draw(self):
        c = self.canv
        def box(x, y, w, h, text, color):
            c.setFillColor(colors.HexColor({BLUE:'#edf4fc', GREEN:'#eef6f1', AMBER:'#fff4e3'}[color]))
            c.setStrokeColor(color)
            c.roundRect(x, y, w, h, 7, fill=1, stroke=1)
            p = Paragraph(text, self.style)
            _, ph = p.wrap(w-16, h-10)
            if ph > h-10:
                raise ValueError('Diagram text exceeds its box')
            p.drawOn(c, x+8, y+(h-ph)/2)
        def arrow(x1,y1,x2,y2):
            import math
            c.setStrokeColor(colors.HexColor('#64758a'))
            c.setLineWidth(1)
            c.line(x1,y1,x2,y2)
            a=math.atan2(y2-y1,x2-x1)
            for d in (-0.45,0.45):
                c.line(x2,y2,x2-6*math.cos(a+d),y2-6*math.sin(a+d))
        ru = self.lang=='RU'
        if self.number==1:
            bw=(W-44)/3; gap=22
            labels = ([
                '<b>Оптимизация</b><br/>Меняем окно CDRH3 трастузумаба',
                '<b>FACS + SPR</b><br/>10 предварительных вариантов',
                '<b>De novo</b><br/>RFD4 создаёт 64 пары VH/VL',
                '<b>Локальный фильтр</b><br/>43 пары проходят',
                '<b>Общий пул: 53</b><br/>Abguard: 52 проходят<br/>Ранжированный черновик: 30'] if ru else [
                '<b>Optimization</b><br/>Change the trastuzumab CDRH3 window',
                '<b>FACS + SPR</b><br/>10 provisional variants',
                '<b>De novo</b><br/>RFD4 generates 64 VH/VL pairs',
                '<b>Local screen</b><br/>43 pairs pass',
                '<b>Combined pool: 53</b><br/>Abguard: 52 pass<br/>Ranked review draft: 30'])
            box(0,120,bw,62,labels[0],BLUE);box(bw+gap,120,bw,62,labels[1],BLUE)
            box(0,24,bw,62,labels[2],GREEN);box(bw+gap,24,bw,62,labels[3],GREEN)
            box(2*(bw+gap),61,bw,84,labels[4],AMBER)
            arrow(bw,151,bw+gap-3,151);arrow(bw,55,bw+gap-3,55)
            arrow(2*bw+gap,151,2*(bw+gap)-3,119);arrow(2*bw+gap,55,2*(bw+gap)-3,87)
        else:
            bw=(W-24)/3
            labels=([
                '<b>Abguard: формат и новизна</b> - допуск, не измеренное связывание',
                '<b>FACS / SPR</b><br/>Проверяем прогноз на измеренных данных',
                '<b>Структура</b><br/>WT: 1/3 и 0/3<br/>Порог 2/3 не достигнут',
                '<b>Термостабильность</b><br/>Tm кандидатов неизвестна',
                '<b>Лаборатория:</b> SPR, конкуренция с трастузумабом, Tm'] if ru else [
                '<b>Abguard: format and novelty</b> - eligibility, not measured binding',
                '<b>FACS / SPR</b><br/>Test predictions against measured data',
                '<b>Structure</b><br/>WT: 1/3 and 0/3<br/>Required 2/3 not reached',
                '<b>Thermal stability</b><br/>Candidate Tm remains unknown',
                '<b>Laboratory:</b> SPR, trastuzumab competition, Tm'])
            box(0,155,W,36,labels[0],BLUE)
            for i,col in enumerate((BLUE,AMBER,AMBER)):
                x=i*(bw+12);box(x,64,bw,70,labels[i+1],col)
                arrow(x+bw/2,155,x+bw/2,137);arrow(x+bw/2,64,x+bw/2,44)
            box(0,5,W,36,labels[4],GREEN)


def build(lang):
    source=ROOT/f'docs/PIPELINE_EXPLAINED_{lang}.md'
    output=ROOT/f'output/pdf/her2-pipeline-explained-{lang.lower()}.pdf'
    output.parent.mkdir(parents=True,exist_ok=True)
    styles={
        'body':ParagraphStyle('body',fontName='Guide',fontSize=10.5,leading=15,textColor=INK,spaceAfter=8),
        'title':ParagraphStyle('title',fontName='GuideBold',fontSize=23,leading=28,textColor=INK,spaceAfter=16),
        'heading':ParagraphStyle('heading',fontName='GuideBold',fontSize=15,leading=20,textColor=BLUE,spaceBefore=14,spaceAfter=9,keepWithNext=True),
        'small':ParagraphStyle('small',fontName='Guide',fontSize=9,leading=12,textColor=INK,spaceAfter=5),
        'diagram':ParagraphStyle('diagram',fontName='Guide',fontSize=9.5,leading=12,textColor=INK),
        'bullet':ParagraphStyle('bullet',fontName='Guide',fontSize=10.5,leading=15,textColor=INK,leftIndent=12,firstLineIndent=-8,spaceAfter=6),
    }
    story=[]; lines=source.read_text().splitlines();i=0;diagram=0
    while i<len(lines):
        line=lines[i].strip()
        if not line: i+=1;continue
        if line.startswith('```mermaid'):
            i+=1
            while i<len(lines) and not lines[i].startswith('```'):i+=1
            diagram+=1;story.extend([Spacer(1,4),Diagram(diagram,lang,styles['diagram']),Spacer(1,8)]);i+=1;continue
        if line.startswith('# '):story.append(Paragraph(inline(line[2:],source),styles['title']));i+=1;continue
        if line.startswith('## '):
            if 'Короткий текст' in line or 'Short speaking' in line:
                story.append(PageBreak())
            story.append(Paragraph(inline(line[3:],source),styles['heading']));i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[-:]+',x) for x in cells):
                    rows.append([Paragraph(inline(x,source),styles['small']) for x in cells])
                i+=1
            table=Table(rows,colWidths=[W*.25,W*.28,W*.47],repeatRows=1,hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#edf4fc')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,0),.7,BLUE),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#d9e1ea'))]))
            story.extend([table,Spacer(1,12)]);continue
        if line.startswith(('**Важно про допуск', '**Eligibility and tracks', 'В FACS модель', 'The FACS model')):
            story.append(PageBreak())
        bullet=line.startswith('- ')
        story.append(Paragraph(('• '+inline(line[2:],source)) if bullet else inline(line,source),styles['bullet' if bullet else 'body']))
        i+=1
    def footer(c,doc):
        c.setFont('Guide',8);c.setFillColor(colors.HexColor('#63758b'))
        c.drawString(48,27,'HER2 • 14 November 2026 • '+lang)
        c.drawRightString(A4[0]-48,27,str(doc.page))
        c.setStrokeColor(colors.HexColor('#d9e1ea'));c.line(48,40,A4[0]-48,40)
    SimpleDocTemplate(str(output),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=42,bottomMargin=54,title=('HER2: пояснение пайплайна' if lang=='RU' else 'HER2: pipeline explained'),author='HER2 hackathon team').build(story,onFirstPage=footer,onLaterPages=footer)
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font-dir',type=Path,default=Path('/System/Library/Fonts/Supplemental'))
    args=parser.parse_args()
    pdfmetrics.registerFont(TTFont('Guide',str(args.font_dir/'Arial.ttf')))
    pdfmetrics.registerFont(TTFont('GuideBold',str(args.font_dir/'Arial Bold.ttf')))
    pdfmetrics.registerFontFamily('Guide',normal='Guide',bold='GuideBold',italic='Guide',boldItalic='GuideBold')
    for lang in ('RU','EN'):print(build(lang))
