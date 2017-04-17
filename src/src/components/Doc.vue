<template>
  <el-row class="doc">
    <el-col :span="6" class="sidebar">
      <div class="doc-nav">
        <h2>Documentation</h2>
        <p><a href="https://github.com/hust220/nsp">Github</a></p>
        <ul>
          <li class="doc-theme" v-for="theme in doc_nav">
            <strong><a :href="theme.href" :class="[theme.href == '#' + $route.path ? 'active' : '']" v-text="theme.text"></a></strong>
            <ul>
              <li class="doc-item" v-for="item in theme.items">
                <a :href="item.href" :class="[item.href == '#' + $route.path ? 'active' : '']" v-text="item.text"></a>
              </li>
            </ul>
          </li>
        </ul>
      </div>
    </el-col>
    <el-col :span="18" :offset="6" class="main">
      <div class="doc-content" v-html="content" :style="'padding-bottom: ' + padding_bottom + 'px'"></div>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'
  import { elementPosition, ScrollToControl } from '../scroll.js'
  import docCache from '../data.js'
  import marked from 'marked'

  export default {
    data() {
      return {
        doc_nav: '',
        content: '',
        padding_bottom: 0
      }
    },

    methods: {
      scroll(item) {
        if (this.imgLoaded()) {
          ScrollToControl(item)
        } else {
          window.setTimeout(() => { this.scroll(item) }, 50)
        }
      },

      fetch_doc_nav() {
        var url = 'static/docs/nav.json'
        this.$http.get(url).then(response => {
          this.doc_nav = response.body
        }, response => {
          //
        })
      },

      set_content(content, f) {
        this.content = content
        this.$nextTick(f)
      },

      fetch_doc_content() {
        var params = this.$route.params
        var theme = (params.theme ? params.theme : 'install')
        var item = (params.item ? params.item : '')
        console.log(docCache)
        console.log(theme + ' ' + item)
        if (docCache.theme && theme === docCache.theme) {
          console.log(1)
          this.scroll(item)
        } else if (theme in docCache) {
          console.log(2)
          this.set_content(docCache[theme], () => {
            this.scroll(item)
            docCache.theme = theme
          })
        } else {
          console.log(3)
          var url = 'static/docs/' + theme + '.md'
          this.$http.get(url).then(response => {
            this.set_content(marked(response.body), () => {
              console.log(item)
              this.scroll(item)
              docCache[theme] = this.content
              docCache.theme = theme
            })
          }, response => {
            //
          })
        }
      },

      imgLoaded() {
        var ls = document.querySelectorAll('.doc-content img')
        var len = ls.length
        for (var i = 0; i < len; i++) {
          var e = ls[i]
          if (!e.complete) return false
        }
        return true
      },

      aboveTopVisualArea(e) {
        try {
          var hw = window.scrollY
          var r = e.href.match(/#\/doc\/.+\/(.+)$/)
          var item = r[1]
          var el = document.getElementById(item)
          var p = elementPosition(el).y
          return p <= hw + 60
        } catch (e) {
          return false
        }
      },

      routeTheme() {
        var params = this.$route.params
        var theme = (params.theme ? params.theme : 'install')
        return theme
      },

      pathTheme(path) {
        try {
          var r = path.match(/#\/doc\/(.+)\/(.+)$/)
          return r[1]
        } catch (e) {
          return ''
        }
      },

      handleScroll() {
        var ls = document.querySelectorAll('.doc-item a')
        var len = ls.length
        var theme = this.routeTheme()
        for (var i = 0; i < len; i++) {
          var e1 = ls[i]
          var e2 = ls[i + 1]
          if (this.pathTheme(e1.href) === theme && this.aboveTopVisualArea(e1) && !this.aboveTopVisualArea(e2)) {
            e1.setAttribute('class', 'active')
          } else {
            e1.setAttribute('class', '')
          }
        }
      },

      handleResized() {
        var h = window.innerHeight
        this.padding_bottom = h - 150
      }

    },

    created() {
      this.fetch_doc_nav()
      this.fetch_doc_content()
      this.handleResized()
      window.addEventListener('scroll', this.handleScroll)
      window.addEventListener('resize', this.handleResized)
    },

    watch: {
      '$route'(to, from) {
        this.fetch_doc_content()
      }
    },

    destroyed() {
      docCache.theme = ''
    }

  }
</script>

<style>
  .sidebar {
    position: fixed;
    font-size: 0.85em;
    overflow-y:scroll;
    height: 100%;
    padding-left: 50px;
  }

  .sidebar .doc-nav {
    padding-bottom: 100px;
  }

  .sidebar ul {
    margin: 10px 15px;
    padding: 0px;
  }

  .sidebar li {
    list-style-type: none;
    margin: 10px 0px;
  }

  .doc-nav >ul {
    margin: 10px 0px;
  }

  .doc-content {
    margin: 35px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  }

  .doc-content pre {
    background-color: #f7f7f7;
    padding: 16px;
  }

  .doc-content code {
    font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  }

  .doc-content table {
    border-collapse: collapse;
  }

  .doc-content tr:nth-child(2n) {
    background-color: #f8f8f8;
  }

  .doc-content td {
    padding: 6px 13px;
    border: 1px solid #ddd;
  }

  .doc-content img {
    max-width: 80%;
  }

  .sidebar a {
    text-decoration: none;
    color: black;
  }

  .sidebar a.active {
    color: blue;
  }

  .sidebar a:hover {
    color: blue;
  }

</style>
