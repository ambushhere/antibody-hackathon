"""Build the HER2 explainer and its overview using ReportLab, Poppler and Pillow."""

import argparse
import math
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageOps
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
W, H = 1200, 675
M = 48
BG = HexColor("#F5F7FA")
INK = HexColor("#132639")
MUTE = HexColor("#475D71")
LINE = HexColor("#D5DFE7")
BLUE = HexColor("#215987")
DARK_BLUE = HexColor("#153A60")
SKY = HexColor("#A9C9DF")
LIGHT_BLUE = HexColor("#E7F0F8")
GREEN = HexColor("#276856")
LEAF = HexColor("#B6D6C6")
LIGHT_GREEN = HexColor("#E7F2EC")
ORANGE = HexColor("#B64C20")
PEACH = HexColor("#FAEBDD")
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
LINUX_FONTS = Path('/usr/share/fonts/truetype/dejavu')

for name, (filename, fallback) in {
    "text": ("Arial.ttf", 'DejaVuSans.ttf'),
    "bold": ("Arial Bold.ttf", 'DejaVuSans-Bold.ttf'),
    "head": ("Trebuchet MS Bold.ttf", 'DejaVuSans-Bold.ttf'),
    "mono": ("Courier New Bold.ttf", 'DejaVuSansMono-Bold.ttf'),
}.items():
    font_path = FONT_DIR / filename
    if not font_path.exists():
        font_path = LINUX_FONTS / fallback
    pdfmetrics.registerFont(TTFont(name, str(font_path)))


class Deck:
    def __init__(self, path):
        self.c = canvas.Canvas(str(path), pagesize=(W, H), pageCompression=1, invariant=1)
        self.c.setTitle("Что мы делаем на AI × Antibodies")
        self.c.setSubject("Мишень HER2, критерии дизайна, бенчмарк и план экспериментов")
        self.page = 0

    def rect(self, x, y, w, h, fill, radius=0, edge=None, width=1):
        self.c.setFillColor(fill)
        if edge:
            self.c.setStrokeColor(edge)
            self.c.setLineWidth(width)
        if radius:
            self.c.roundRect(x, H-y-h, w, h, radius, fill=1, stroke=bool(edge))
        else:
            self.c.rect(x, H-y-h, w, h, fill=1, stroke=bool(edge))

    def line(self, x1, y1, x2, y2, color=LINE, width=1, dash=None):
        self.c.setStrokeColor(color)
        self.c.setLineWidth(width)
        self.c.setLineCap(1)
        self.c.setDash(dash or [])
        self.c.line(x1, H-y1, x2, H-y2)
        self.c.setDash([])

    def arrow(self, x1, y1, x2, y2, color=MUTE, width=1.7, tip=7):
        self.line(x1, y1, x2, y2, color, width)
        angle = math.atan2(y2 - y1, x2 - x1)
        for offset in (-0.56, 0.56):
            x = x2 - tip * math.cos(angle + offset)
            y = y2 - tip * math.sin(angle + offset)
            self.line(x, y, x2, y2, color, width)

    def path(self, start, commands, color, width=2):
        self.c.setStrokeColor(color)
        self.c.setLineWidth(width)
        self.c.setLineCap(1)
        self.c.setLineJoin(1)
        p = self.c.beginPath()
        p.moveTo(start[0], H-start[1])
        for cmd in commands:
            if len(cmd) == 2:
                p.lineTo(cmd[0], H-cmd[1])
            else:
                x1, y1, x2, y2, x3, y3 = cmd
                p.curveTo(x1, H-y1, x2, H-y2, x3, H-y3)
        self.c.drawPath(p)

    def dot(self, x, y, radius, color, edge=None, width=1):
        self.c.setFillColor(color)
        if edge:
            self.c.setStrokeColor(edge)
            self.c.setLineWidth(width)
        self.c.circle(x, H-y, radius, fill=1, stroke=bool(edge))

    def tx(self, x, y, content, size=20, color=INK, font="text", align="left"):
        """Draw a line with y approximately at the top of its capital letters."""
        span = pdfmetrics.stringWidth(content, font, size)
        left = {"left": x, "center": x-span/2, "right": x-span}[align]
        if left < 0 or left+span > W+.05 or y < 0 or y+size > H:
            raise ValueError(f"Page {self.page}: text outside page: {content}")
        self.c.setFont(font, size)
        self.c.setFillColor(color)
        self.c.drawString(left, H-y-size*.77, content)
        return span

    def body(self, x, y, content, width, size=21, color=MUTE, font="text", leading=None):
        leading = leading or size*1.32
        for paragraph in content.split("\n"):
            line = ""
            for word in paragraph.split():
                trial = f"{line} {word}" if line else word
                if line and pdfmetrics.stringWidth(trial, font, size) > width:
                    self.tx(x, y, line, size, color, font)
                    y += leading
                    line = word
                else:
                    line = trial
                if pdfmetrics.stringWidth(line, font, size) > width:
                    raise ValueError(f"Page {self.page}: word too wide: {word}")
            self.tx(x, y, line, size, color, font)
            y += leading
        return y

    def link(self, x, y, label, url, size=16):
        w = self.tx(x, y, label, size, BLUE)
        self.c.linkURL(url, (x, H-y-size-2, x+w, H-y+2), relative=0, thickness=0)

    def start(self, section, title=None, sub=None):
        self.page += 1
        self.rect(0, 0, W, H, BG)
        self.rect(M, 37, 22, 4, ORANGE, 2)
        self.tx(M+33, 31, section, 16, BLUE, "bold")
        self.tx(W-M, 31, "Осло  /  14 ноября 2026", 15, MUTE, align="right")
        if title:
            self.tx(M, 70, title, 46, INK, "head")
        if sub:
            self.tx(M, 131, sub, 21, MUTE)

    def end(self, source):
        self.line(M, 622, W-M, 622)
        self.tx(M, 640, source, 14, MUTE)
        self.tx(W-M, 637, f"{self.page:02d} / 10", 16, BLUE, "bold", "right")
        self.c.showPage()

    def note(self, x, y, w, text, color=BLUE, fill=LIGHT_BLUE, size=19):
        self.rect(x, y, w, 46, fill, 10)
        self.body(x+14, y+15, text, w-28, size, color, "bold")


def receptor(d, x, y):
    d.rect(x+27, y+324, 450, 35, LIGHT_GREEN, 10)
    for dx in range(39, 469, 16):
        d.dot(x+dx, y+330, 2.5, LEAF)
        d.dot(x+dx, y+354, 2.5, LEAF)
    d.line(x+314, y+114, x+278, y+351, GREEN, 8)
    d.rect(x+264, y+321, 28, 44, GREEN, 11)
    d.rect(x+242, y+368, 72, 42, LEAF, 17)
    for i, (dx, dy, w, h) in enumerate([
        (277, 57, 91, 63), (260, 132, 72, 52),
        (246, 197, 81, 55), (240, 264, 77, 48),
    ]):
        d.rect(x+dx, y+dy, w, h, GREEN if i == 3 else LEAF, 17)
        d.tx(x+dx+w/2, y+dy+h/2-8, ("I", "II", "III", "IV")[i], 19,
             white if i == 3 else GREEN, "bold", "center")
    d.rect(x+238, y+275, 10, 25, HexColor("#80B79B"), 4)
    d.c.saveState()
    d.c.translate(x+120, H-y-246)
    d.c.rotate(63)
    d.c.translate(-x-120, -(H-y-246))
    igg(d, x+120, y+177, .74)
    d.c.restoreState()
    d.path((x+193,y+281), [(x+207,y+271,x+211,y+296,x+233,y+287)], ORANGE, 6)
    d.tx(x+26, y+29, "Схема контакта", 15, MUTE, "bold")
    d.tx(x+411, y+76, "HER2", 24, GREEN, "head", "center")
    d.line(x+369, y+92, x+382, y+92, GREEN)
    d.tx(x+52, y+167, "Трастузумаб", 22, BLUE, "bold")
    d.line(x+119,y+192,x+126,y+219,BLUE)
    d.tx(x+342, y+279, "Эпитоп", 20, GREEN, "bold")
    d.path((x+251,y+287), [(x+331,y+287)], GREEN, 1.2)
    d.tx(x+36, y+334, "Мембрана", 16, GREEN, "bold")
    d.tx(x+340, y+379, "Внутри клетки", 16, GREEN)


def cover(d):
    d.start("AI × Antibodies")
    d.tx(M, 124, "Что мы делаем", 67, INK, "head")
    d.tx(M, 202, "на AI × Antibodies", 60, INK, "head")
    d.body(M, 300, "Проектируем варианты антитела, которые должны узнавать тот же участок HER2, что и трастузумаб.", 495, 28)
    d.body(M, 432, "HER2 проходит через мембрану клетки. Трастузумаб узнаёт участок его наружного домена IV, рядом с мембраной.", 477, 22)
    d.note(M, 542, 474, "Исследовательские эксперименты без подачи")
    d.rect(625, 106, 527, 482, white, 20, LINE)
    receptor(d, 634, 118)
    d.tx(661, 559, "Условная схема, без атомных координат", 16, MUTE)
    d.end("Мишень и референс [1, 2, 6]")


def igg(d, x, y, scale=1):
    def p(dx, dy):
        return x+dx*scale, y+dy*scale

    def rod(x1, y1, x2, y2, color, thickness=22):
        d.line(*p(x1, y1), *p(x2, y2), color, thickness*scale)

    for sign in (-1, 1):
        rod(sign*13, 86, sign*13, 181, DARK_BLUE)
        rod(sign*13, 87, sign*91, 15, DARK_BLUE)
        rod(sign*100, 29, sign*55, 73, SKY)
        rod(sign*64, 39, sign*91, 15, BLUE)
        rod(sign*100, 29, sign*81, 47, HexColor("#6B9CC0"))
    rod(-11, 88, 11, 88, DARK_BLUE, 7)


def site(d):
    d.start("Мишень и антитело", "Один сайт связывания, шесть петель",
            "У IgG два сайта. В каждом сходятся вариабельные домены тяжёлой (VH) и лёгкой (VL) цепей.")
    d.rect(48, 178, 325, 424, white, 18, LINE)
    d.rect(396, 178, 756, 424, white, 18, LINE)
    d.tx(72, 205, "Целое антитело IgG", 25, INK, "head")
    igg(d, 202, 286, .98)
    d.dot(299, 310, 26, BG, ORANGE, 2)
    d.line(297, 304, 306, 297, BLUE, 9)
    d.line(300, 318, 311, 308, SKY, 9)
    d.line(325, 310, 373, 310, ORANGE, 1.5)
    d.arrow(373, 310, 408, 310, ORANGE, 1.5)
    d.tx(220, 400, "Fc", 18, MUTE, "bold")
    d.dot(82, 524, 5, DARK_BLUE)
    d.tx(96, 515, "2 тяжёлые цепи", 20)
    d.dot(82, 560, 5, SKY)
    d.tx(96, 551, "2 лёгкие цепи", 20)
    d.tx(423, 205, "Один кончик IgG крупным планом", 25, INK, "head")
    d.rect(442, 271, 327, 62, LEAF, 14)
    d.tx(605, 286, "HER2, домен IV", 21, GREEN, "bold", "center")
    d.line(523, 330, 701, 330, GREEN, 7)
    d.line(769, 330, 811, 330, GREEN, 1.5)
    d.tx(827, 298, "Эпитоп", 23, GREEN, "head")
    d.tx(827, 331, "Участок на мишени HER2.", 19, MUTE)
    d.rect(445, 409, 159, 106, LIGHT_BLUE, 14)
    d.rect(610, 409, 159, 106, LIGHT_BLUE, 14)
    d.line(460, 358, 750, 358, BLUE, 1.3, [4, 4])
    for i, (x, y, label) in enumerate([
        (459, 382, "H1"), (505, 373, "H2"), (551, 350, "H3"),
        (625, 365, "L1"), (672, 381, "L2"), (719, 375, "L3"),
    ]):
        color = ORANGE if i == 2 else BLUE
        d.path((x, 416), [(x+5, y, x+27, y, x+34, 416)], color, 5)
        d.tx(x+17, 424, label, 15, color, "bold", "center")
    d.tx(525, 469, "VH", 27, BLUE, "head", "center")
    d.tx(689, 469, "VL", 27, BLUE, "head", "center")
    d.path((752, 358), [(784, 358), (807, 396)], BLUE, 1.4)
    d.tx(827, 381, "Паратоп", 23, BLUE, "head")
    d.body(827, 414, "Контактная поверхность антитела. Её в основном образуют CDR, петли на VH и VL.", 293, 19)
    d.tx(445, 550, "3 петли на VH + 3 на VL", 21, BLUE, "bold")
    d.tx(827, 529, "CDRH3", 23, ORANGE, "head")
    d.tx(827, 564, "Третья петля тяжёлой цепи.", 19, MUTE)
    d.end("Схема доменов и петель [3, 5, 6]")


def domain(d, x, y, name, fill=LIGHT_BLUE, color=BLUE, w=54, h=41):
    d.rect(x, y, w, h, fill, 8)
    d.tx(x+w/2, y+h/2-7, name, 17, color, "bold", "center")


def formats(d):
    d.start("Форматы и треки", "Что именно мы проектируем",
            "Формат определяет, какие домены входят в конструкцию. Трек определяет, что разрешено менять.")
    columns = [
        ("IgG", "Полное антитело с двумя сайтами и областью Fc."),
        ("VHH", "Один вариабельный домен тяжёлой цепи."),
        ("Парные VH / VL", "Две отдельные цепи. В задании это «Fab»."),
        ("scFv", "VH и VL соединены гибким линкером в одну цепь."),
    ]
    for i, (title, description) in enumerate(columns):
        x = 48+280*i
        d.rect(x, 178, 264, 185, white, 14, LINE)
        d.tx(x+16, 196, title, 24, INK, "head")
        d.body(x+16, 303, description, 232, 18, leading=22)
        if i == 0:
            igg(d, x+132, 227, .35)
        elif i == 1:
            domain(d, x+101, 242, "VHH", w=62)
        else:
            domain(d, x+65, 240, "VH")
            domain(d, x+144, 240, "VL", SKY, DARK_BLUE)
            if i == 3:
                d.path((x+107, 277), [(x+120, 307, x+132, 225, x+151, 244)], BLUE, 2.5)
    d.body(48, 379, "Классический Fab включает также константные домены CH1 и CL. Пара VH/VL без них называется Fv.", 1080, 19)
    d.rect(48, 426, 540, 174, PEACH, 15)
    d.tx(69, 447, "Оптимизация трастузумаба", 25, ORANGE, "head")
    d.tx(69, 484, "W G G D G F Y A M D", 28, ORANGE, "mono")
    d.body(69, 531, "Меняем это окно из 10 остатков внутри CDRH3. Остальной VH, вся VL и формат сохраняются.", 493, 20)
    d.rect(608, 426, 544, 174, LIGHT_GREEN, 15)
    d.tx(629, 447, "De novo, тот же эпитоп", 25, GREEN, "head")
    d.body(629, 489, "Ищем достаточно отличающееся антитело к тому же участку HER2. Проверка новизны учитывает формат и идёт по базе Adaptyv.", 475, 20)
    d.end("Форматы и правила треков [1, 2, 3, 4]")


def binding(d):
    d.start("Лабораторные критерии", "Три свойства, три измерения",
            "Сравниваем с трастузумабом в сопоставимом формате и при одинаковых условиях.")
    for x in (48, 423, 798):
        d.rect(x, 180, 354, 377, white, 17, LINE)
    for x, title, subtitle in [
        (48, "Связывание", "Измеряем KD методом SPR"),
        (423, "Эпитоп", "Измеряем конкуренцию"),
        (798, "Термостабильность", "Измеряем температуру Tm"),
    ]:
        d.tx(x+19, 203, title, 28, INK, "head")
        d.tx(x+19, 244, subtitle, 19, MUTE)

    d.rect(68, 282, 133, 58, LIGHT_BLUE, 10)
    d.rect(249, 282, 132, 58, LIGHT_GREEN, 10)
    d.tx(134, 305, "Антитело + HER2", 15, BLUE, "bold", "center")
    d.tx(315, 305, "Комплекс", 17, GREEN, "bold", "center")
    d.arrow(211, 303, 239, 303, BLUE, 1.6, 5)
    d.arrow(239, 319, 211, 319, BLUE, 1.6, 5)
    d.tx(68, 364, "KD = koff / kon", 26, BLUE, "head")
    d.tx(253, 371, "модель 1:1", 16, MUTE)
    d.body(68, 407, "Меньше KD означает сильнее связывание в равновесии. ΔΔG связывания показывает, как изменилась свободная энергия связывания.", 311, 20)

    d.rect(451, 326, 298, 37, LIGHT_GREEN, 10)
    d.line(577, 326, 624, 326, GREEN, 8)
    d.rect(463, 282, 120, 29, LIGHT_BLUE, 8)
    d.tx(523, 289, "Трастузумаб", 15, BLUE, "bold", "center")
    d.rect(614, 282, 123, 29, PEACH, 8)
    d.tx(675, 289, "Кандидат", 15, ORANGE, "bold", "center")
    d.arrow(569, 312, 581, 322, BLUE, 1.6, 5)
    d.arrow(632, 312, 620, 322, ORANGE, 1.6, 5)
    d.tx(600, 338, "HER2", 18, GREEN, "bold", "center")
    d.body(443, 394, "Мешает ли кандидат трастузумабу занять участок HER2? Конкуренция указывает на перекрытие или близость участков, включая взаимное экранирование.", 313, 20)

    d.path((833, 314), [(836,290,865,281,870,297), (880,318,849,328,847,308),
                      (843,290,900,286,893,311), (893,325,862,340,847,332)], BLUE, 4)
    d.arrow(915, 305, 987, 305, MUTE, 1.6, 6)
    d.arrow(987, 319, 915, 319, MUTE, 1.6, 6)
    d.path((1000,309), [(1012,285,1021,337,1035,310), (1048,278,1067,335,1080,312),
                       (1095,290,1100,327,1120,304)], BLUE, 4)
    d.tx(868, 351, "Свёрнуто", 16, BLUE, "bold", "center")
    d.tx(1058, 351, "Развёрнуто", 16, BLUE, "bold", "center")
    d.body(818, 394, "Tm характеризует переход при нагреве. Более высокий Tm означает, что белок теряет свёрнутую форму при более высокой температуре.", 313, 20)
    d.note(48, 568, 1104, "Аффинность оцениваем по KD. ipTM оценивает уверенность в интерфейсе, ESM likelihood описывает последовательность.", size=17)
    d.end("Методы оценки [1, 2, 4]")


def data(d):
    d.start("Экспериментальные данные", "FACS и SPR отвечают на разные вопросы",
            "Метки получены на вариантах трастузумаба. У новых молекул эти свойства ещё предстоит измерить.")
    d.rect(48, 177, 534, 373, white, 18, LINE)
    d.rect(614, 177, 538, 373, white, 18, LINE)
    d.tx(72, 202, "FACS", 32, BLUE, "head")
    d.tx(638, 202, "SPR", 32, GREEN, "head")
    d.tx(72, 249, "Массовый скрининг", 23, INK, "bold")
    d.tx(638, 249, "Количественная аффинность", 23, INK, "bold")
    for i, (label, color, amount) in enumerate([
        ("low", LIGHT_BLUE, 1), ("mid", SKY, 2), ("high", BLUE, 3)
    ]):
        x = 75+157*i
        d.rect(x, 301, 143, 47, color, 9)
        for j in range(amount):
            d.dot(x+17+j*10, 325, 3, white if i == 2 else BLUE)
        d.tx(x+125, 316, label, 18, white if i == 2 else BLUE, "bold", "right")
    for i, (value, caption) in enumerate([
        ("768", "строк анализа"), ("190", "строк с KD"), ("158", "вариантов в оценке"),
    ]):
        x = 639+166*i
        d.tx(x, 296, value, 41, GREEN, "head")
        d.tx(x, 348, caption, 16, MUTE)
        if i < 2:
            d.arrow(x+105, 319, x+138, 319, GREEN)
    d.body(72, 382, "Сортировка клеток по флуоресцентному сигналу даёт три класса. Для ML задача бинарная: выделить high среди mid и low.", 475, 22)
    d.body(638, 382, "Датчик измеряет связывание и диссоциацию во времени. KD количественно описывает аффинность. Среди 190 строк есть повторные измерения исходного варианта WT.", 477, 22)
    d.note(48, 568, 1104, "Модель обучается на FACS. SPR проверяет, переносится ли порядок её оценок на аффинность.", size=19)
    d.end("Данные организаторов и отбор SPR [2, 3, 4]")


def benchmark(d):
    d.start("Проверка ML", "Что проверили на исторических данных",
            "524 346 FACS-вариантов. Использованы опубликованные сплиты remove_overlap по клонотипам.")
    for i, (label, count, detail) in enumerate([
        ("Train", "367 044", "Обучаем веса логистической модели."),
        ("Validation", "78 651", "Выбираем модель по AP и фиксируем выбор."),
        ("Test / holdout", "78 651", "Измеряем качество после фиксации выбора."),
    ]):
        x = 48+375*i
        d.rect(x, 179, 354, 142, LIGHT_BLUE, 13)
        d.tx(x+15, 198, label, 19, BLUE, "bold")
        d.tx(x+15, 225, count, 35, INK, "head")
        d.body(x+15, 274, detail, 322, 17, leading=20)
        if i < 2:
            d.arrow(x+358, 249, x+370, 249, MUTE, tip=5)
    d.tx(48, 338, "Между сплитами 0 совпадений CDRH3 и центров клонотипов (групп похожих петель).", 20, MUTE)
    d.rect(48, 378, 540, 222, white, 15, LINE)
    d.rect(608, 378, 544, 222, white, 15, LINE)
    d.tx(68, 396, "FACS test  /  high против остальных", 22, BLUE, "bold")
    d.tx(68, 439, "0,9859", 45, INK, "head")
    d.tx(303, 439, "0,9737", 45, INK, "head")
    d.tx(68, 491, "ROC AUC", 17, BLUE, "bold")
    d.tx(303, 491, "Average precision", 17, BLUE, "bold")
    d.body(68, 534, "Логистическая регрессия с признаками отдельных позиций и их пар.", 496, 19)
    d.tx(629, 396, "SPR  /  158 подходящих CDRH3", 22, GREEN, "bold")
    d.tx(629, 439, "0,35", 43, INK, "head")
    d.tx(866, 439, "0,61", 43, INK, "head")
    d.tx(629, 491, "FACS → SPR", 17, GREEN, "bold")
    d.tx(866, 491, "Обучение на SPR", 17, GREEN, "bold")
    d.body(629, 529, "Ранговая корреляция. Для SPR-модели близкие петли удержаны в одной группе при nested CV. Новых измерений нет.", 503, 18)
    d.end("Метрики проверены в материалах проекта [4]. Новые молекулы в этот бенчмарк не входят.")


def gates(d):
    d.start("Допуск кандидатов", "Почему локальный фильтр строже",
            "Источники задают разные границы для оптимизации. Фильтр выбирает более строгую трактовку.")
    d.rect(48, 178, 1104, 164, white, 14, LINE)
    for x, caption in [(67, "Источник"), (523, "Идентичность CDRH3 к WT"), (905, "До выданных данных")]:
        d.tx(x, 195, caption, 18, MUTE, "bold")
    d.line(67, 222, 1132, 222)
    for y, label, identity, distance in [
        (236, "Слайды Apheris, слайд 6", "≤ 70%", "≥ 2 мутаций"),
        (278, "Текст задания, стр. 3 и 7", "< 70%", "> 2 мутаций"),
    ]:
        d.tx(67, y, label, 22, INK, "bold")
        d.tx(523, y, identity, 24, BLUE, "head")
        d.tx(905, y, distance, 22, BLUE, "bold")
    d.tx(67, 317, "На стр. 2 текста также указано ≥ 2. WT означает исходный трастузумаб.", 15, MUTE)
    d.rect(48, 359, 540, 164, PEACH, 15)
    d.rect(608, 359, 544, 164, LIGHT_BLUE, 15)
    d.tx(69, 379, "≥ 4 замены к WT в окне", 27, ORANGE, "head")
    d.tx(629, 379, "≥ 3 правки до ближайшей петли", 26, BLUE, "head")
    d.body(69, 420, "3 замены оставляют 7 из 10 совпадений (70%). 4 оставляют 6 (60%). Идентичность CDRH3 отдельно проверена во всех пяти схемах immunum.*", 497, 19)
    d.body(629, 420, "Локально считаем расстояние Левенштейна до выданных CDRH3. Официальная база Adaptyv и итоговая классификация недоступны.", 499, 19)
    d.body(48, 538, "Общие условия: формат VHH, scFv или VH/VL, ≤ 250 остатков, без непарных цистеинов, immunum ≥ 0,1.", 1104, 19, INK)
    d.body(48, 576, "* IMGT, Kabat, Chothia, Martin, AHo. Проверены 504 записи выборки для локального отбора.", 1104, 16)
    d.end("Правила и локальная трактовка [1, 2, 4, 5]")


def experiments(d):
    d.start("Текущие эксперименты", "От генерации к отбору кандидатов",
            "Статистика предлагает последовательности. RFD4 предлагает последовательности и геометрию.")
    for i, (number, title, caption) in enumerate([
        ("10 000", "Статистический пул", "Варианты окна CDRH3 по частотам и мутациям."),
        ("4", "RFD4 для оптимизации", "Варианты трастузумаба с изменённым окном."),
        ("4", "Редизайн CDR", "Петли на известном каркасе. Новизна de novo пока открыта."),
    ]):
        x = 48+375*i
        d.rect(x, 178, 354, 174, white, 15, LINE)
        d.tx(x+17, 195, number, 45, ORANGE, "head")
        d.tx(x+17, 252, title, 22, INK, "bold")
        d.body(x+17, 290, caption, 316, 18, leading=22)
    d.tx(48, 370, "Ещё 8 дизайнов всей последовательности, две настройки. Новизна по базе Adaptyv не проверена.", 20, MUTE)
    d.tx(48, 415, "Порядок вычислительных проверок", 26, INK, "head")
    for i, (title, caption, fill) in enumerate([
        ("Фильтр", "Формат, нумерация, новизна.", LIGHT_BLUE),
        ("FACS и SPR", "Два сигнала для вариантов трастузумаба.", LIGHT_BLUE),
        ("Кофолдинг", "Модель комплекса, контакты с HER2, контроль WT.", LIGHT_GREEN),
        ("Отбор", "Разнообразие и риски. Проверенного предиктора Tm нет.", PEACH),
    ]):
        x = 48+281*i
        d.rect(x, 459, 260, 95, fill, 12)
        d.tx(x+14, 474, title, 20, INK, "bold")
        d.body(x+14, 507, caption, 231, 17, leading=19)
        if i < 3:
            d.arrow(x+265, 507, x+276, 507, MUTE, tip=5)
    d.body(48, 573, "Высокий скор помогает выбрать, что проверять. Связывание, конкуренция и термостабильность требуют измерений.", 1104, 19, INK)
    d.end("Артефакты кампании и протокол отбора [4]")


def lab(d):
    d.start("Перспективная проверка", "Как проверить гипотезу отбора",
            "При малом числе анализов важен выбор кандидатов с разными основаниями для успеха.")
    d.rect(48, 177, 416, 423, DARK_BLUE, 18)
    d.tx(73, 211, "4", 96, white, "head")
    d.body(76, 323, "гарантированных лабораторных слота на команду, общих для двух треков", 346, 24, white, "bold")
    d.line(75, 440, 437, 440, SKY, 1)
    d.tx(76, 463, "До 20 в резервном списке", 22, white, "bold")
    d.body(76, 503, "Так устроен хакатон. Наша работа ограничена экспериментами, конкурсной подачи нет.", 343, 20, white)
    d.rect(489, 177, 663, 423, white, 18, LINE)
    steps = [
        ("Зафиксировать прогноз", "Записать кандидатов, происхождение, данные, критерии успеха и неопределённости до измерений."),
        ("Измерить вместе с референсом", "Синтезировать выбранные варианты. Измерить KD, конкуренцию с трастузумабом и Tm."),
        ("Сопоставить прогнозы с результатом", "Проверить, обогащает ли отбор пул удачными вариантами. Четырёх измерений мало для точной оценки эффекта."),
    ]
    for i, (title, caption) in enumerate(steps):
        y = 206+111*i
        d.dot(525, y+11, 14, LIGHT_BLUE)
        d.tx(525, y+4, str(i+1), 17, BLUE, "bold", "center")
        d.tx(553, y, title, 23, INK, "head")
        d.body(553, y+36, caption, 564, 19)
    d.body(514, 542, "Гипотеза: прозрачный отбор по нескольким сигналам повышает шанс найти удачный вариант.", 599, 19, BLUE, "bold")
    d.end("Лабораторные слоты и резерв [1, 2]. План исследований [4].")


def sources(d):
    d.start("Источники", "Откуда взяты правила, цифры и схемы")
    d.rect(48, 142, 1104, 460, white, 18, LINE)
    entries = [
        ("[1]", "2026-11_Oslo_hackathon_Apheris_challenge.pdf",
         "Слайды 5–8. Треки, критерии допуска, распределение слотов и экспериментальная оценка."),
        ("[2]", "ai-antibodies-hackathon-apheris.pdf",
         "Стр. 2–4, 7–8. Форматы, конфликт порогов, данные, слоты и лабораторные методы."),
        ("[3]", "Документация организаторов, разделы 17-her2-data и 03-antibodies",
         "Данные HER2 и устройство антител. В этом снимке доступны названия разделов."),
        ("[4]", "docs/BENCHMARK.md  ·  docs/CAMPAIGN.md",
         "Сплиты, метрики FACS/SPR и правила отбора. Общие результаты в results/experiment-20261114."),
        ("[5]", "immunum, документация по нумерации и границам CDR",
         "pypi.org/project/immunum/"),
        ("[6]", "Cho et al., Nature (2003), PDB 1N8Z",
         "rcsb.org/structure/1N8Z"),
    ]
    for i, (n, title, caption) in enumerate(entries):
        y = 163+72*i
        d.tx(67, y+1, n, 19, BLUE, "bold")
        d.tx(117, y, title, 22, INK, "bold")
        if i == 4:
            d.link(117, y+31, caption, "https://pypi.org/project/immunum/", 18)
        elif i == 5:
            d.tx(117, y+31, "Структура HER2 с Fab трастузумаба. Иллюстрации в презентации условные.", 18, MUTE)
            d.link(780, y+31, caption, "https://www.rcsb.org/structure/1N8Z", 18)
        else:
            d.tx(117, y+31, caption, 18, MUTE)
        if i < 5:
            d.line(117, y+59, 1130, y+59)
    d.end("AI × Antibodies, Осло")


def make_overview(pages, output):
    tile_w, tile_h, gap = 360, 203, 12
    columns = 2
    rows = math.ceil(len(pages)/columns)
    sheet = Image.new("RGB", (columns*tile_w+(columns+1)*gap,
                              rows*tile_h+(rows+1)*gap), "#DAE3ED")
    for i, image in enumerate(pages):
        with Image.open(image) as page:
            thumb = ImageOps.contain(page.convert("RGB"), (tile_w, tile_h), Image.Resampling.LANCZOS)
        x = gap+(i % columns)*(tile_w+gap)
        y = gap+(i // columns)*(tile_h+gap)
        sheet.paste(thumb, (x, y))
    sheet.save(output, optimize=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-dir", type=Path, help="Retain per-slide PNGs in this directory.")
    args = parser.parse_args()
    out = ROOT / "output" / "explainer"
    out.mkdir(parents=True, exist_ok=True)
    pdf = out / "her2-intuition.pdf"
    d = Deck(pdf)
    for slide in (cover, site, formats, binding, data, benchmark, gates, experiments, lab, sources):
        slide(d)
    d.c.save()
    with tempfile.TemporaryDirectory(prefix="her2-explainer-") as temp:
        render_dir = args.render_dir or Path(temp)
        render_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(["pdftoppm", "-r", "96", "-png", str(pdf), str(render_dir / "page")], check=True)
        pages = [render_dir / f"page-{i:02}.png" for i in range(1, d.page+1)]
        make_overview(pages, out / "overview.png")


if __name__ == "__main__":
    main()
