<template>
  <el-row class="api">
    <el-col :span="6" class="sidebar">
      <div class="api-nav">
        <h2><a href="#/api">API</a></h2>
        <p><a href="https://github.com/hust220/nsp">Github</a></p>
        <ul>
          <li class="api-theme" v-for="theme in api_nav">
            <strong><a :href="theme.href" :class="[theme.href == '#' + $route.path ? 'active' : '']" v-text="theme.text"></a></strong>
            <ul>
              <li class="api-item" v-for="item in theme.items">
                <a :href="item.href" :class="[item.href == '#' + $route.path ? 'active' : '']" v-text="item.text"></a>
              </li>
            </ul>
          </li>
        </ul>
      </div>
    </el-col>
    <el-col :span="18" :offset="6" class="main">
      <div class="api-content" v-html="content" :style="'padding-bottom: ' + padding_bottom + 'px'"></div>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'
  import { elementPosition, ScrollToControl } from '../scroll.js'
  import { docCache, apiCache } from '../data.js'
  import marked from 'marked'

  export default {
    data() {
      return {
        api_nav: '',
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

      fetch_api_nav() {
        var url = 'static/apis/nav.json'
        this.$http.get(url).then(response => {
          this.api_nav = response.body
        }, response => {
          //
        })
      },

      set_content(content, f) {
        this.content = content
        this.$nextTick(f)
      },

      fetch_api_content() {
        // var params = this.$route.params
        var theme = this.routeTheme()
        var item = this.routeItem()
        // var theme = (params.theme ? params.theme : 'install')
        // var item = (params.item ? params.item : '')
        console.log(apiCache)
        console.log(theme + ' ' + item)
        if (apiCache.theme && theme === apiCache.theme) {
          console.log(1)
          this.scroll(item)
        } else if (theme in apiCache) {
          console.log(2)
          this.set_content(apiCache[theme], () => {
            this.scroll(item)
            apiCache.theme = theme
          })
        } else {
          console.log(3)
          var url = 'static/apis/' + theme + '.md'
          this.$http.get(url).then(response => {
            this.set_content(marked(response.body), () => {
              console.log(item)
              this.scroll(item)
              apiCache[theme] = this.content
              apiCache.theme = theme
            })
          }, response => {
            //
          })
        }
      },

      imgLoaded() {
        var ls = document.querySelectorAll('.api-content img')
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
          var r = e.href.match(/#\/api\/.+\/(.+)$/)
          var item = r[1]
          var el = document.getElementById(item)
          var p = elementPosition(el).y
          return p <= hw + 60
        } catch (e) {
          return false
        }
      },

      routeItem() {
        var params = this.$route.params
        var item = (params.item ? params.item : '')
        return item
      },

      routeTheme() {
        var params = this.$route.params
        var theme = (params.theme ? params.theme : 'index')
        return theme
      },

      pathTheme(path) {
        try {
          var r = path.match(/#\/api\/(.+)\/(.+)$/)
          return r[1]
        } catch (e) {
          return ''
        }
      },

      handleScroll() {
        var ls = document.querySelectorAll('.api-item a')
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
      this.fetch_api_nav()
      this.fetch_api_content()
      this.handleResized()
      window.addEventListener('scroll', this.handleScroll)
      window.addEventListener('resize', this.handleResized)
    },

    watch: {
      '$route'(to, from) {
        this.fetch_api_content()
      }
    },

    destroyed() {
      apiCache.theme = ''
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

  .sidebar .api-nav {
    padding-bottom: 100px;
  }

  .sidebar ul {
    margin: 10px 15px;
    padding: 0px;
  }

  .sidebar .api-nav >ul {
    margin-left: 0px;
  }

  .sidebar li {
    list-style-type: none;
    margin: 10px 0px;
  }

  .api-content {
    margin: 35px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  }

  .api-content pre {
    background-color: #f7f7f7;
    padding: 16px;
  }

  .api-content code {
    font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  }

  .api-content table {
    border-collapse: collapse;
  }

  .api-content tr:nth-child(2n) {
    background-color: #f8f8f8;
  }

  .api-content td {
    padding: 6px 13px;
    border: 1px solid #ddd;
  }

  .api-content img {
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