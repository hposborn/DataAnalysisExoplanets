
import re

def preprocess(md: str) -> str:

# 1. Handle line-by-line [+] fragments safely, preserving indentation & lists
    new_lines = []
    for line in md.split('\n'):
        if '[+]' in line:
            # Capture leading whitespace (indentation)
            indent_len = len(line) - len(line.lstrip())
            indent = line[:indent_len]
            stripped = line.strip()

            # Remove [+] from the stripped line content
            content_without_plus = stripped.replace('[+]', '').strip()

            # Check if it's a bullet list item (e.g., starts with * or -)
            if content_without_plus.startswith(('* ', '- ', '+ ')):
                bullet = content_without_plus[:2]
                text = content_without_plus[2:].strip()
                # Rebuild preserving indentation, bullet, and wrapping text in span
                clean_line = f'{indent}{bullet}<span class="fragment">{text}</span>'
            else:
                # Regular paragraph or text line
                clean_line = f'{indent}<span class="fragment">{content_without_plus}</span>'

            new_lines.append(clean_line)
        else:
            new_lines.append(line)
    md = '\n'.join(new_lines)# columns: <!-- c -->, <!-- | -->, <!-- c. -->
    md = re.sub(r"<!--\s*c\s*-->", r'\n<div class="cols" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*ct\s*-->", r'\n<div class="cols cols-top-align" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*cb\s*-->", r'\n<div class="cols cols-bottom-align" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*\|\s*-->", r'\n</div><div markdown="1">\n', md)
    md = re.sub(r'<!--\s*c.\s*-->', r'\n</div></div>\n', md)
    md = re.sub(r"<!--\s*c21\s*-->", r'\n<div class="cols21" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*c12\s*-->", r'\n<div class="cols12" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*c31\s*-->", r'\n<div class="cols31" markdown="1"><div markdown="1">\n', md)
    md = re.sub(r"<!--\s*c111\s*-->", r'\n<div class="cols111" markdown="1"><div markdown="1">\n', md)

# single close: </div> (for boxes, footnotes, etc.)
    md = re.sub(r"<!--\s*\.\s*-->", r'\n</div>\n', md)

# width: <!-- w100% -->, <!-- w50% -->
    md = re.sub(r"<!--\s*w([0-9]+)%\s*-->", r'<!-- .element style="width:\1%" -->', md)
    md = re.sub(r"<!--\s*w([0-9]+)\s*-->", r'<!-- .element style="width:\1%" -->', md)

# absolute position with width: <!-- absolute 100 200 10 -->
    md = re.sub(r"<!--\s*absolute ([0-9]+) ([0-9]+) ([0-9]+)\s*-->", r'<!-- .element style="position: absolute; transform: translate(-50%, -50%); top:\1%; left:\2%; width:\3%" -->', md)

# vertical space: <!-- vspace2.5 -->
    md = re.sub(r"<!--\s*vspace([0-9.]+)\s*-->", r'\n<p style="margin-top: \1em;"></p>\n', md)

# left/right float <!-- float-left -->, <!-- float-right -->
    md = re.sub(r"<!--\s*float-left\s*-->", r'<!-- .element style="float: left" -->', md)
    md = re.sub(r"<!--\s*float-right\s*-->", r'<!-- .element style="float: right" -->', md)

# box: <!-- box -->, <!-- . -->
    md = re.sub(r"<!--\s*box\s*-->", r'\n<div class="box" markdown="1">\n', md)

# footnote: <!-- footnote -->, <!-- . -->
    md = re.sub(r"<!--\s*footnote\s*-->", r'\n<div class="footnote" markdown="1">\n', md)

# math: <!-- e -->, <!-- e. -->
    md = re.sub(r"<!--\s*e\s*-->", r'\n<div markdown="1">\\[\\begin{aligned}', md)
    md = re.sub(r"<!--\s*e.\s*-->", r'\\end{aligned}\\]</div>', md)

# fragments: <!-- f0 -->, <!-- f1 -->
    md = re.sub(r"<!--\s*f([0-9]+)\s*-->", r'<!-- .element: class="fragment" data-fragment-index="\1" -->', md)

# fit-text: <!-- fit -->
    md = re.sub(r"<!--\s*fit\s*-->", r'<!-- .element class="r-fit-text" -->', md)

# auto-animate: <!-- anim0 -->, <!-- anim -->
    md = re.sub(r"<!--\s*anim\s*-->", r'<!-- .element data-auto-animate -->', md)
    md = re.sub(r"<!--\s*anim0\s*-->", r'<!-- .element data-auto-animate data-auto-animate-restart -->', md)

# left/right text-align: <!-- text-left -->, <!-- text-right -->
    md = re.sub(r"<!--\s*text-left\s*-->", r'<!-- .element style="text-align: left" -->', md)
    md = re.sub(r"<!--\s*text-right\s*-->", r'<!-- .element style="text-align: right" -->', md)

#  Handle [q] by wrapping the text in a custom styled div or box.
    # We use a CSS class like 'question-box' so you can style it in your CSS file.
    # If someone writes [+]\n[q] text, we can also handle it, or use a specific syntax like [q+]

    # Simple block replacement for [q]
    # We convert [q] followed by text into a styled div
    md = re.sub(
        r'\[q\]\s*(.*)',
        r'<div class="question-box">\1</div>',
        md
    )
    # Helper for local vs external media (automatically handles images OR videos)
    def create_media_tag(match, extra_attrs):
        path = match.group(1).strip().strip('"\'')
        if path.startswith('http'):
            src = path
        else:
            src = f'../slide_data/slide_images/{path}'

        # Check if the file extension is a video format
        video_extensions = ('.mp4', '.webm', '.ogg', '.mov')
        if path.lower().endswith(video_extensions):
            # Render a centered HTML5 video tag with controls, autoplay, loop, and mute (safe for browser autoplay policies)
            return f'<div style="text-align: center; width: 100%;"><video src="{src}" controls autoplay loop muted {extra_attrs}></video></div>'
        else:
            # Render normal image tag
            return f'<div style="text-align: center; width: 100%;"><img src="{src}" {extra_attrs} /></div>'

    # 4. Handle [i=...] macro (safe default sizing, prevents overflowing the page)
    md = re.sub(
        r'\[i=[\'"]?([^\]\'"]+)[\'"]?\]',
        lambda m: create_media_tag(m, 'style="max-height: 450px; max-width: 100%; object-fit: contain; display: block; margin: 0 auto;"'),
        md
    )

    # 5. Handle [ism=...] image macro (350px height, with or without quotes)
    md = re.sub(
        r'\[ism=[\'"]?([^\]\'"]+)[\'"]?\]',
        lambda m: create_media_tag(m, 'style="height: 250px;"'),
        md
    )

    # 6. Handle [ifull=...] image macro (full screen / max slide height coverage)
    md = re.sub(
        r'\[ifull=[\'"]?([^\]\'"]+)[\'"]?\]',
        lambda m: create_media_tag(m, 'style="width: 100%; max-height: 620px; object-fit: contain; display: block; margin: 0 auto;"'),
        md
    )

    # 6. Automatic Auto-Numbered Contents & Parts Generator
    part_pattern = re.compile(r'\[P\s+["\']([^"\']+)["\']\]')

    # Find all titles first to build the [contents] list
    titles = part_pattern.findall(md)

    if '[contents]' in md:
        if titles:
            contents_html = '<ul class="contents-list">\n'
            for idx, title in enumerate(titles, start=1):
                contents_html += f'  <li><a href="#/part-{idx}"><strong>Part {idx}:</strong> {title}</a></li>\n'
            contents_html += '</ul>'
        else:
            contents_html = '<p style="color: gray;"><em>(No [P...] parts found in this file)</em></p>'

        md = md.replace('[contents]', contents_html)

    # Replace each [P "title"] tag with an auto-incremented number, heading, and anchor ID
    part_counter = [0]
    def replace_part(match):
        part_counter[0] += 1
        num = part_counter[0]
        title = match.group(1)
        return f'<h2 id="part-{num}" class="slide-part-title">Part {num}: {title}</h2>'

    md = part_pattern.sub(replace_part, md)
    return md

#
# import re
#
# def preprocess(md: str) -> str:
#
#
# # Short way of making the fragment appear:
#     md = md.replace("[+]", "<!-- .element: class=\"fragment\" -->\n")
#
# # columns: <!-- c -->, <!-- | -->, <!-- c. -->
#     md = re.sub(r"<!--\s*c\s*-->", r'\n<div class="cols" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*ct\s*-->", r'\n<div class="cols cols-top-align" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*cb\s*-->", r'\n<div class="cols cols-bottom-align" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*\|\s*-->", r'\n</div><div markdown="1">\n', md)
#     md = re.sub(r'<!--\s*c.\s*-->', r'\n</div></div>\n', md)
#     md = re.sub(r"<!--\s*c21\s*-->", r'\n<div class="cols21" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*c12\s*-->", r'\n<div class="cols12" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*c31\s*-->", r'\n<div class="cols31" markdown="1"><div markdown="1">\n', md)
#     md = re.sub(r"<!--\s*c111\s*-->", r'\n<div class="cols111" markdown="1"><div markdown="1">\n', md)
#
# # single close: </div> (for boxes, footnotes, etc.)
#     md = re.sub(r"<!--\s*\.\s*-->", r'\n</div>\n', md)
#
# # width: <!-- w100% -->, <!-- w50% -->
#     md = re.sub(r"<!--\s*w([0-9]+)%\s*-->", r'<!-- .element style="width:\1%" -->', md)
#     md = re.sub(r"<!--\s*w([0-9]+)\s*-->", r'<!-- .element style="width:\1%" -->', md)
#
# # absolute position with width: <!-- absolute 100 200 10 -->
#     md = re.sub(r"<!--\s*absolute ([0-9]+) ([0-9]+) ([0-9]+)\s*-->", r'<!-- .element style="position: absolute; transform: translate(-50%, -50%); top:\1%; left:\2%; width:\3%" -->', md)
#
# # vertical space: <!-- vspace2.5 -->
#     md = re.sub(r"<!--\s*vspace([0-9.]+)\s*-->", r'\n<p style="margin-top: \1em;"></p>\n', md)
#
# # left/right float <!-- float-left -->, <!-- float-right -->
#     md = re.sub(r"<!--\s*float-left\s*-->", r'<!-- .element style="float: left" -->', md)
#     md = re.sub(r"<!--\s*float-right\s*-->", r'<!-- .element style="float: right" -->', md)
#
# # box: <!-- box -->, <!-- . -->
#     md = re.sub(r"<!--\s*box\s*-->", r'\n<div class="box" markdown="1">\n', md)
#
# # footnote: <!-- footnote -->, <!-- . -->
#     md = re.sub(r"<!--\s*footnote\s*-->", r'\n<div class="footnote" markdown="1">\n', md)
#
# # math: <!-- e -->, <!-- e. -->
#     md = re.sub(r"<!--\s*e\s*-->", r'\n<div markdown="1">\\[\\begin{aligned}', md)
#     md = re.sub(r"<!--\s*e.\s*-->", r'\\end{aligned}\\]</div>', md)
#
# # fragments: <!-- f0 -->, <!-- f1 -->
#     md = re.sub(r"<!--\s*f([0-9]+)\s*-->", r'<!-- .element: class="fragment" data-fragment-index="\1" -->', md)
#
# # fit-text: <!-- fit -->
#     md = re.sub(r"<!--\s*fit\s*-->", r'<!-- .element class="r-fit-text" -->', md)
#
# # auto-animate: <!-- anim0 -->, <!-- anim -->
#     md = re.sub(r"<!--\s*anim\s*-->", r'<!-- .element data-auto-animate -->', md)
#     md = re.sub(r"<!--\s*anim0\s*-->", r'<!-- .element data-auto-animate data-auto-animate-restart -->', md)
#
# # left/right text-align: <!-- text-left -->, <!-- text-right -->
#     md = re.sub(r"<!--\s*text-left\s*-->", r'<!-- .element style="text-align: left" -->', md)
#     md = re.sub(r"<!--\s*text-right\s*-->", r'<!-- .element style="text-align: right" -->', md)
#
# #  Handle [q] by wrapping the text in a custom styled div or box.
#     # We use a CSS class like 'question-box' so you can style it in your CSS file.
#     # If someone writes [+]\n[q] text, we can also handle it, or use a specific syntax like [q+]
#
#     # Simple block replacement for [q]
#     # We convert [q] followed by text into a styled div
#     md = re.sub(
#         r'\[q\]\s*(.*)',
#         r'<div class="question-box">\1</div>',
#         md
#     )
#     return md
