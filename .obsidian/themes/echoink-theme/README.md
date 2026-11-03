# echoink-theme

An offline Obsidian companion theme for EchoInk. Download the separately provided theme from [echoink.cn/theme](https://echoink.cn/theme). Version 0.1.0, Obsidian 1.11.4 or newer.

Copy this folder to your vault’s `.obsidian/themes/`, then choose **echoink-theme** in Settings → Appearance → Themes. It is not listed in the community directory yet.

The theme follows Obsidian’s light/dark mode. EchoInk’s Appearance setting selects Minimal (`minimal`), Highlighter (`highlighter`), or Manifesto (`manifesto`). Without EchoInk, with no selection, or with an unknown identifier, Minimal is used. Disabling EchoInk restores Minimal. Minimal is free; Highlighter and Manifesto require active EchoInk PRO. Your preference is retained while inactive and restored when PRO becomes available again.

With Obsidian’s **Readable line length** enabled, reading and editing content use 80% of the note pane, capped at 1200px and centered, with at least about 16px of side space in narrow panes. Properties and embedded content share that width. Turning the option off keeps Obsidian’s native full-width behavior. EchoInk’s own pages are unaffected.

All fonts and icons are embedded. No network requests are made. Embedded font licenses and source credits are in `OFL-LICENSES.txt`; property and task icons use IconPark SVGs (Apache-2.0, see `ICON-LICENSE.txt`).

Reading view alone uses tilted paper elements. Live Preview and source mode remain upright. Ordinary image embeds always remain separate. Use an explicit gallery callout for image groups (one image, two side by side, or three with a large image on the left). Add `|row` to the callout type for three equal columns. Live Preview rendering and reading view share the layout; revealing the Markdown source temporarily uses the native source layout. With EchoInk enabled, empty groups, groups containing other content, and groups with more than three images show a three-line error instead of their contents. Click the error in Live Preview to reveal the source and correct it. Without EchoInk, unsupported groups retain ordinary visible content.

Native content folding is retained. Theme-owned transitions last 180ms and are disabled with reduced motion. Custom task characters change appearance, while native checkbox interactions remain unchanged.

```md
> [!gallery]
> ![[one.png]]
> ![[two.png]]
> ![[three.png]]

> [!gallery|row]
> ![[one.png]]
> ![[two.png]]
> ![[three.png]]
```

With EchoInk and the gallery-capable theme active, desktop Live Preview provides image selection circles and a grouping toolbar. All platforms provide Insert image gallery, Group selected images, and Split image gallery commands. Select 2–3 consecutive image-only lines (blank lines are allowed). Original paths, dimensions, and aliases are preserved; each action can be undone once. No circles are added on phones or tablets.

Inside a gallery, width hints such as `|508` do not constrain a tile: each photo fills its allocated cell and is cropped to fit. The original Markdown dimensions and aliases stay intact, so splitting the group restores the original size hints.


Write one image per quoted line (recommended); space-separated or adjacent embeds on the same line and blank `>` lines also work. Counts select the layout: one full-width image (280px high), two equal columns (280px), three with one large left image and two stacked right images (360px total), or `|row` for three equal columns (220px). `row` with two images still means two columns.

To set the height, use `[!gallery|300]` or `[!gallery|row|300]`; `[!gallery|300|row]` also works. Heights are clamped to 120–800px. Photos fill and crop to their cells; per-image resize handles are hidden inside a group. Live Preview adds only a dashed frame and a layout/height label, keeping image geometry identical to reading view. The source opens normally when the cursor enters the group.

Theme-only installation supports default layouts and heights for canonical image-only Markdown callouts. Manual heights and format errors require EchoInk; theme-only fallback never hides malformed content.
